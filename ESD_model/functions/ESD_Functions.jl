function adapt_global_coupling(hparams,N::Int64,W::Matrix{Float64},lags::Matrix{Float64},h,t::Float64,u::Vector{Float64})
    @inbounds for ii = 1:N
        @inbounds for jj = 1:N
            if W[jj,ii]  > 0.0
                if lags[jj,ii] == 0.0
                    W[ii,jj]+= 0.001*u[ii]*(u[jj] - h(hparams,t-1.0;idxs=jj))
                else
                    W[ii,jj]+= 0.001*h(hparams,t-lags[jj,ii];idxs=ii)*(u[jj] - h(hparams,t-1.0;idxs=jj))
                end


                if W[ii,jj] > SC[ii,jj] + SC[ii,jj]*0.05
                    W[ii,jj] = SC[ii,jj] + SC[ii,jj]*0.05
                elseif W[ii,jj] < SC[ii,jj] - SC[ii,jj]*0.05
                    W[ii,jj] = SC[ii,jj] - SC[ii,jj]*0.05
                end

        
                #W[ii,jj] += W[ii,jj]*(1-exp(-abs(SC[ii,jj] - W[ii,jj])/0.01))*dw
            end
        end
     

    end

    @inbounds for k=1:N #Maintaining symmetry in the weights between regions
        @inbounds for l = k:N
                 W[k,l] = (W[k,l]+W[l,k])/2
                 W[l,k] = W[k,l]
        end
    end

    W .= W./maximum(W)
    
    return W
end

function f(x::Float64,θ::Float64,a::Float64,θj::Float64,aj::Float64) 
    return (1. /(1+exp(-(aj/(1+a))*(x-(θj + θ))))) -  (1. /(1+exp((aj*θj)/(1+a))))
end

function k(a::Float64,θj::Float64,aj::Float64) 
    return ((1+exp((aj*θj)/(1+a))) /(1+exp((aj*θj)/(1+a))))
end





function adapt_local_func(w1,w2,w3,w4,w5,w6,w7,wsum,E,S,D,h,hparams,t,i,N,c)


    hE = E - h(hparams,t-1.0;idxs = i)
    hS = S - h(hparams,t-1.0,idxs = i+N)
    hD = D - h(hparams,t-1,idxs=i+2*N)
    
    w1 += c*E*hE
    w2 += c*S*hE
    w3 += c*D*hE
    w4 += c*E*hS
    w6 += c*E*hD

    wsum_new = w1+w2+w3+w4+w5+w6+w7

    w1 = (wsum*w1)/(wsum_new)
    w2 = (wsum*w2)/(wsum_new)
    w3 = (wsum*w3)/(wsum_new)
    w4 = (wsum*w4)/(wsum_new)
    w6 = (wsum*w6)/(wsum_new)


    return w1,w2,w3,w4,w6
end


	

function getModelFC(BOLD_TRIALS,nTrials,tstart,step,size)
    modelFC = []
    for i = 1:nTrials
        if i == 1
            modelFC = get_FC(BOLD_TRIALS[:,:,i],tstart,step,size)/nTrials
        else
            modelFC += get_FC(BOLD_TRIALS[:,:,i],tstart,step,size)/nTrials
        end
    end
    return modelFC
end

function remove_zeros_W!(W,N)
    indx_zeros = 0 
    break_cond = false
    for i = 1:size(W,3)
        if W[:,:,i] == zeros(N,N)
            indx_zeros = i
            break
        end
    end

    return W[:,:,1:indx_zeros-1]
end
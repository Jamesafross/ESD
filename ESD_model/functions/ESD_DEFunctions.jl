function ESD_DE(du,u,h,p,t)
    hparams = p

    @unpack modelparams,weights,network,stimparams,helpers,plasticity_opts,solver_opts = ESD
 
    @unpack w1,w2,w3,w4,w5,w6,w7,PE,PI,θE,θI,aE,aI,η,τE,τI,τISP,ρ = modelparams
    @unpack W,lags,N = network
    @unpack w1v,w2v,w3v,w4v,w6v,w_sum = weights
    @unpack plast_time,learning_rate = plasticity_opts

    @unpack stimWindow,stimNodes,stimStr,Tstim = stimparams
    @unpack wmat,d,current_window,adapt_time,tprev= helpers


    make_hist_mat2_threads!(h,W,u,hparams,N,lags,t,wmat)
    
    sum!(d,wmat)

    @inbounds Threads.@threads for i = 1:N

        E = u[i]
        S = u[i+N]
        D = u[i+2N]

        if ((t == round(adapt_time,digits=3) && (t != tprev)) || (t - tprev > 0.01)) && t >= plast_time 
            ESD.weights.w1v[i],ESD.weights.w2v[i],
            ESD.weights.w3v[i], ESD.weights.w4v[i], 
            ESD.weights.w6v[i] = adapt_local_func(w1v[i],w2v[i],w3v[i],w4v[i],w5,w6v[i],w7,w_sum,E,S,D,h,hparams,t,i,N,learning_rate)
            if i == N
                ESD.helpers.adapt_time += 0.01
                ESD.helpers.tprev = maximum([tprev,t])
                ESD.helpers.adapt_counter += 1

                if mod(ESD.helpers.adapt_counter,100)==0
                    ESD.weights_save.w1[ESD.weights_save.count,:] = ESD.weights.w1v
                    ESD.weights_save.w2[ESD.weights_save.count,:] = ESD.weights.w2v
                    ESD.weights_save.w3[ESD.weights_save.count,:] = ESD.weights.w3v
                    ESD.weights_save.w4[ESD.weights_save.count,:] = ESD.weights.w4v
                    ESD.weights_save.w6[ESD.weights_save.count,:] = ESD.weights.w6v
                    ESD.weights_save.count += 1
                end
                
            end

           

        end


        s = stim(t,i,stimparams,current_window)


        du[i] = (1/τE)*(-E + (k(w3v[i]*D,θE,aE)-E)*f(w1v[i]*E + η*d[i] + s + PE, w2v[i]*S,w3v[i]*D,θE,aE))
        du[i+N] =(1/τI)*( -S + (k(0.,θI,aI)-S)*f(w4v[i]*E + PI, 0., 0.,θI,aI))
        du[i+2N] =(1/τI)*( -D + (k(0.,θI,aI)-D)*f(w5*E+PI, w6v[i]*S + w7*D, 0.0 ,θI,aI))
     

    end
end


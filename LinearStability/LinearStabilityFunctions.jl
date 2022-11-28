function f(x,β,θ)
    return 1/(1+exp(-β*(x-θ)))
end

function Df(x,β,θ)
    return (β*(exp(-β*(x-θ))))/((exp(-β*(x-θ)) + 1)^2)
end

function μFunc(p,λ,v)
    μ = 0.

    for i = 1:N; for j = 1:N;
        μ += SC[i,j]*exp(-λ*dist[i,j])*v[:,p][i]*v[:,p][j]
    end;end;

    return μ
end

function AFunc(params,p,Ebar,Ibar,λ,v)
    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,β,η,σ = params
    τ_mat = [τE 0. ; 0. τI]
    I2 = [1. 0. ; 0. 1.]
    wloc = [cEE cEI; cIE cII]
    μ = μFunc(p,λ,v)

    Df1 = Df(cEE*Ebar[p] + cEI*Ibar[p] +  sum(SC[p,:].*Ebar[1:N]) ,β,θE)
    Df2 = Df(cIE*Ebar[p] + cII*Ibar[p],β,θI)
    return -τ_mat*(I2 - [Df1 ; Df2].*(wloc .+ [μ 0. ; 0. 0. ]))

end



function network_steady_state!(F,x,params)
    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,β,η,σ = params
    for i = 1:N
        F[i] = -x[i] + f(cEE*x[i] + cEI*x[i+N] + sum(SC[i,:].*x[1:N]) ,β,θE)
        F[i+N] = -x[i+N] + f(cIE*x[i] + cII*x[i+N],β,θI)
    end
end
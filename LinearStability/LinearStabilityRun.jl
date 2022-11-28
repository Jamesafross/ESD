using NLsolve,DifferentialEquations,Random,StatsBase,Parameters,JLD,LinearAlgebra,Plots
include("LinearStabilityFunctions.jl")
include("LinearStabilityParams.jl")
include("NetworkSetup.jl")


c = 13.
const SC,dist,lags,N = networksetup(c;digits=3,nSC=1,nFC=1,N=140,normalise=false)

eigs = eigen(SC)
λ = eigs.values
v = eigs.vectors
v = v./maximum(v)



WCp = WCparams()


f!(F,x) = network_steady_state!(F,x,WCp)


steadyState  = nlsolve(f!,0.1.*ones(2N)).zero
while  ( sum(steadyState .< 0) > 0 )  
    global steadyState  = nlsolve(f!,0.5.*rand(2N)).zero
end


Ebar = steadyState[1:N]
Ibar = steadyState[N+1:2N]
p=1

AFunc(WCp,p,Ebar,Ibar,λ[p],v)

function λ_p1!(F,x,p,params,Ebar,Ibar,v)
    A = AFunc(params,p,Ebar,Ibar,x[1]+im*x[2],v)
    F[1] = real(det(A))
    F[2] = imag(det(A))
end

function λ_n1!(F,x,p,params,Ebar,Ibar,v)
    A = AFunc(params,p,Ebar,Ibar,x[1]+im*x[2],v)
    F[1] = real(det(A))
    F[2] = imag(det(A))
end


λ_p = zeros(N)

λ_n = zeros(N)

for i = 1:N
    println(i)
    λ_p!(F,x) = λ_p1!(F,x,i,WCp,Ebar,Ibar,v)
    λ_n!(F,x) = λ_n1!(F,x,i,WCp,Ebar,Ibar,v)

    flag = false
    while flag == false 
        test = nlsolve(λ_p!,[randn(),randn()])
        flag = test.f_converged
    
        λ_p[i] = test.zero[1]
    end

   
end
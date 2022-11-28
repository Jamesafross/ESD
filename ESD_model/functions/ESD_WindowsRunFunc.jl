function ESD_windows_run()
    @unpack modelparams,weights,balloonparams,network,stimparams,helpers,solver_opts = ESD
 
    @unpack w1,w2,w3,w4,w5,w6,w7,PE,PI,θE,θI,aE,aI,η,τE,τI,τISP,ρ = modelparams
    @unpack W,lags,N = network
   

    @unpack stimWindow,stimNodes,stimStr,Tstim = stimparams
    @unpack wmat,d,current_window = helpers

    @unpack tWindow, nWindow = solver_opts
    BOLD_saveat = collect(0:1.0:tWindow)
    size_out = length(BOLD_saveat)
    BOLD_out = zeros(N,size_out,nWindow)
    
    for j = 1:nWindow

        if j >= ESD.plasticity_opts.plast_start_win
            ESD.plasticity_opts.plast_time = 0.0
            ESD.helpers.adapt_time = 0.0
        end

    
        ESD.helpers.current_window = j
       
        println("working on window . . . ",j)
        if j == 1
            ESD.IC = zeros(3N)
            hparams = ESD.IC

            h = h1 
        else
            ESD.IC = sol[:,end]
            iStart = findfirst(sol.t .> tWindow - 1.2)
            u_hist = make_uhist(sol.t[iStart:end] .- sol.t[end],sol[:,iStart:end])
            hparams = u_hist
            h = h2
        end

        tspan = (0.0,tWindow)
        
        p = hparams
        
    
        prob = DDEProblem(ESD_DE, ESD.IC, h, tspan, p)
        println("solving...")
        adp_times = collect(0.0:0.01:tWindow)
        global sol = solve(prob,MethodOfSteps(BS3()),maxiters = 1e20,saveat=0.01,tstops=adp_times)

        BalloonIn= make_In(sol.t,sol[1:N,:])
        tspanB = (sol.t[1],sol.t[end])
        balloonParams = balloonparams,BalloonIn,N
        if j == 1
            b0 =  cat(zeros(N),ones(3N),dims=1)
        else
            b0 = endBM
        end
        println(". . . Running Balloon Model")

        global out,endBM = runBalloon(b0,balloonParams,tspanB,BOLD_saveat,N)
        
        if j == 1
            ESD.bold_out = out
        else
            ESD.bold_out = cat(ESD.bold_out,out[:,2:end],dims=2)
    
        end

    end

end



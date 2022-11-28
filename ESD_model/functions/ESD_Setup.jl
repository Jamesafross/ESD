function setup()
   #load data and make struct & dist matrices
   
   println("setting up Wilson Cowan Model")
   ESDp = ESDparams()
   c = 7000
   constant_delay = 0.005

   global SC,dist,lags,N,FC,missingROIs = networksetup(c,constant_delay;N=140,adj=false)
   
   wmat = zeros(N,N)
   d=zeros(N)
   @unpack w1,w2,w3,w4,w5,w6,w7 = ESDp
   network_params = networkParameters(SC,dist,round.(lags,digits=5),N)
   pweights = plastic_weights(ones(N).*w1,w2*ones(N),w3*ones(N),w4*ones(N),w6*ones(N),w1+w2+w3+w4+w5+w6+w7)
   p_opts = PlasticityOptions()
   p_opts.plast_time = 100.0
   p_opts.plast_start_win = 2
   weights_save = weightsSave(zeros(N),zeros(N),zeros(N),zeros(N),zeros(N),1)

   ESD = ESDstruct(ESDp,
                  pweights,
                  network_params,
                  balloonModelParameters(),
                  zeros(4N),
                  StimParams(),
                  SolverOptions(),
                  Helpers(wmat,d,1,p_opts.plast_time,p_opts.plast_time,0),
                  p_opts,
                  zeros(N),
                  weights_save
   )
  
   return ESD
end



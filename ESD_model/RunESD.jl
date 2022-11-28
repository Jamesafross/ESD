include("functions/ESD_InitSetup.jl")
using RollingFunctions

global ESD = setup()

ESD.modelparams.PE=1.65
ESD.modelparams.η=0.1
ESD.solver_opts.tWindow = 60
ESD.solver_opts.nWindow = 300
ESD.plasticity_opts.plast_start_win = 30
ESD.stimparams.stimWindow = 150
ESD.plasticity_opts.learning_rate = 0.01
ESD.stimparams.stimNodes = [39,13,36,19,15,35]
ESD.stimparams.stimStr = -1.0
ESD.stimparams.stim= true

nsave = (ESD.solver_opts.nWindow-28)*ESD.solver_opts.tWindow
ESD.weights_save = weightsSave(zeros(nsave,N),zeros(nsave,N),zeros(nsave,N),zeros(nsave,N),zeros(nsave,N),1)


@time ESD_windows_run()

modelFC = getModelFC(ESD.bold_out,1,100,100,300) 
#controlFC,missingROIs = get_mean_all_functional_data(;ROI=140,type="control")
#FC_fit_to_data_mean = zeros(size(modelFC,3))
#FC_STIM = load("$PROGDIR/_data/FC_stimulated_times.jld","FC_stimulated_times")
#for j = 1:size(modelFC,3)
 #   FC_fit_to_data_mean[j] = fit_r(modelFC[:,:,j],controlFC)
#end

#plot(FC_fit_to_data_mean)


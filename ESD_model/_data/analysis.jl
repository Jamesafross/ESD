using JLD,Plots,RollingFunctions,Graphs,SimpleWeightedGraphs,StatsBase,Parameters
WORKDIR = "/Users/james/ESD/ESD_model"
include("$WORKDIR/../GlobalFunctions/_includes.jl")
include("_functions.jl")
nn0 = 2
nn = 3
stim1 = load("$WORKDIR/_data/weights_STIM_ROI140_$(nn).jld","weights_STIM_ROI140_$(nn)")
control = load("$WORKDIR/_data/weights_NO_STIM_ROI140_$nn.jld","weights_NO_STIM_ROI140_$nn")
SC,dist,lags,N,FC,missingROIs = networksetup(13000,0.05;N=140)




D = SimpleWeightedGraph(dist)
S = SimpleWeightedGraph( 1 ./ SC)
sp_S = johnson_shortest_paths(S).dists
sp_D = johnson_shortest_paths(D).dists
sp_D = sp_D[4,:]
sp_S = sp_S[4,:]

di=100*(abs.((stim1 .- control))./control)
diff_weights = mean(di,dims=3)[:,:]'
percent_diff =diff_weights

percent_diff[isnan.(percent_diff)] .= 0 


rollmean_pd = roll_mean_mat(percent_diff[:,7000:end],1)'
rollmean_pd[39,:] .= 0

rollmean_pos = rollmean_pd[rollmean_pd .> 0]
m = mean(rollmean_pos)
st = std(rollmean_pos)



Y = zeros(size(rollmean_pd))
Y .= rollmean_pd

Y[Y .< m + 1.5st] .= 0


ROIs,times = get_start_times(Y)

times_f = zeros(length(times))
dists = zeros(length(times))
strs = zeros(length(times))

global counter = 1
for i in ROIs
    times_f[counter] = times[counter]
    dists[counter] = sp_D[ROIs[counter]]
    strs[counter] = sp_S[ROIs[counter]]
    global counter += 1
end

println("correlation (distance) = ", corspearman(dists,times_f))
println("correlation (stregnths) = ", corspearman(strs,times_f))

heatmap(Y,color=:jet)



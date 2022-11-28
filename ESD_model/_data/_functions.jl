function get_start_times(y)
    ROIs = Array{Int64}[]
    times = Array{Float64}[]

    for i = 1:size(y,1)
        time = findfirst(x ->x > 0, y[i,:])
        if (findfirst(x ->x > 0, y[i,:]) == nothing) == false
            ROIs = cat(ROIs,i,dims=1)
            times = cat(times,time,dims=1)
        end
    end
    return ROIs, times
end


function threshold_!(Y,x)

    for i=1:size(Y,1)
        m = mean(Y[i,:])
        st = std(Y[i,:])
        V = Y[i,:]
        V[V .< m + x*st] .= 0
        Y[i,:] = V
    end


end
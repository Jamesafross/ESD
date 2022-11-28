@with_kw mutable struct ESDparams{R}
    w1::R = 15.
    w2::R  = 12.
    w3::R  = 15.
    w4::R  = 14.
    w5::R = 20.
    w6::R = 16.
    w7::R= 8.
    θE::R  = 4.
    θI::R  = 3.7
    aE::R = 1.3
    aI::R = 2.0
    η::R  = 0.12
    σ::R  = 0.005
    τE::R  = 0.01
    τI::R  =0.02
    τx::R = .01
    PE::R =2.
    PI::R = 1.
    τISP::R=0.25
    ρ::R=0.30
end


mutable struct networkParameters
    W::Matrix{Float64}
    dist::Matrix{Float64}
    lags::Matrix{Float64}
    N::Int64
end

mutable struct plastic_weights
    w1v::Array{Float64}
    w2v::Array{Float64}
    w3v::Array{Float64}
    w4v::Array{Float64}
    w6v::Array{Float64}
    w_sum::Float64
end

mutable struct Helpers
    wmat::Array{Float64}
    d::Array{Float64}
    current_window::Real
    adapt_time::Float64
    tprev::Float64
    adapt_counter::Int
end

@with_kw mutable struct PlasticityOptions
    plast_start_win::Real = 200
    plast_time::Real = 0.01
    learning_rate::Real=0.005
end


@with_kw mutable struct StimParams
    stim::Bool = false
    stimWindow::Real = 20
    stimNodes::Vector{Real} = [39]
    stimStr::Real = -1.
    Tstim::Vector{Real} = [10,40]
end

@with_kw mutable struct SolverOptions
    tWindow::Real = 6000.
    nWindow::Real = 1
end


@with_kw mutable struct weightsSave
    w1::Array{Real}
    w2::Array{Real}
    w3::Array{Real}
    w4::Array{Real}
    w6::Array{Real}
    count::Int = 1
end

mutable struct ESDstruct
    modelparams::ESDparams
    weights::plastic_weights
    network::networkParameters
    balloonparams::balloonModelParameters
    IC::Vector{Float64}
    stimparams::StimParams
    solver_opts::SolverOptions
    helpers::Helpers
    plasticity_opts::PlasticityOptions
    bold_out
    weights_save::weightsSave
end




	

using NPZ
using CSV
using Plots
using FFTW
using Optim
using Statistics
using DataFrames
using Distributed
using OffsetArrays
using SparseArrays
using StaticArrays
using SharedArrays
using ProgressMeter
using LinearAlgebra
using BenchmarkTools
using DelimitedFiles

# include custom code
include("./Test/runtests.jl")
include("./Source/ActivePolymer.jl")
using .ActivePolymer
using .ActivePolymer.CorrelationMatrices

profile = npzread("../Share/ABidentities_blobel2021_chr2_35Mb_60Mb.npy") .==1

ratio   = 5.974;
diff    = (ratio - 1)/(ratio+1);
Tₐ      = 1 .+ diff * 2*(profile .- 0.5);
T₀      = 1 .+ 0.0  * 2*(profile .- 0.5);

Nₜ = 1000;
τs  = exp.(range(-6,18,Nₜ));

ΔR²reference_mean = fill(0.0,Nₜ)
ΔR²reference_n500 = fill(0.0,Nₜ)
@showprogress for (id, τ) in enumerate(τs)
    msd = ActivePolymer.Transform.Forward.compute_mean_squared_traveled_distance(
        T₀ |> diagm, 
        ActivePolymer.Jacobian.Discrete.J₀,
        τ,
        fourier_type=ActivePolymer.Methods.FastFourier.DCT
    )
    ΔR²reference_mean[id] = msd |> mean
    ΔR²reference_n500[id] = msd[500]
end

ΔR²mean = fill(0.0,Nₜ)
ΔR²hot_mean  = fill(0.0,Nₜ)
ΔR²cold_mean = fill(0.0,Nₜ)
ΔR²hot_n350  = fill(0.0,Nₜ)
ΔR²cold_n500 = fill(0.0,Nₜ)
@showprogress for (id, τ) in enumerate(τs)
    msd = ActivePolymer.Transform.Forward.compute_mean_squared_traveled_distance(
        Tₐ |> diagm, 
        ActivePolymer.Jacobian.Discrete.J₀,
        τ,
        fourier_type=ActivePolymer.Methods.FastFourier.DCT
    );
    ΔR²mean[id] = msd |> mean
    ΔR²hot_mean[id]  = msd[profile] |> mean
    ΔR²cold_mean[id] = msd[.!profile] |> mean
    ΔR²hot_n350[id]  = msd[350]
    ΔR²cold_n500[id] = msd[500]    
end

data = DataFrame(
    τs = τs,
    reference_mean = ΔR²reference_mean,
    reference_n500 = ΔR²reference_n500,
    mean = ΔR²mean,
    hot_mean = ΔR²hot_mean,
    cold_mean = ΔR²cold_mean,
    hot_n350 = ΔR²hot_n350,
    cold_n500 = ΔR²cold_n500
);
CSV.write("discrete_diffusion_msd.csv", data)
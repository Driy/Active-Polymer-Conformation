using NPZ
using CSV
using FFTW
using HDF5
using Optim
using PyPlot
using PyCall
using DataFrames
using Statistics
using Distributed
using OffsetArrays
using SparseArrays
using StaticArrays
using SharedArrays
using LinearAlgebra
using ProgressMeter
using BenchmarkTools
using DelimitedFiles
using ImageFiltering

#
cmasher = pyimport("cmasher")

# include custom code
include("./Test/runtests.jl")
include("./Source/ActivePolymer.jl")
using .ActivePolymer
using .ActivePolymer.CorrelationMatrices

association = [
    "Deq1" 1
    "bcomps_2x" 2
    "bcomps_3x" 3
    "bcomps_4x" 4
    "bcomps_5x" 5
    "bcomps_7x" 7
    "bcomps_10x" 10
    "bcomps_19x" 19
    "bcomps_26x" 26
    "bcomps_39x" 39
];
names = association[:,1];
values = association[:,2];
mappedvalues = 2(values .- 1) ./ (values .+1);

analysis_type = "optimization"
folder = ["output_", analysis_type] |> join
mkpath(folder)

# setup reference
R, ΔR = ActivePolymer.Optimization.Interface.load_data(names[1]);

# Three-parameter fits
modeltype  = ActivePolymer.Optimization.Model.Full
jacmodule  = ActivePolymer.Jacobian.Discrete
n          = 3

parameters = ActivePolymer.Optimization.Interface.fit_mechanics(
    ΔR, modeltype=modeltype, jacmodule=jacmodule, n=n, padding=0.85)

# export Jacobian
jacobian   = jacmodule.J(parameters.minimizer[2:end]..., n)
    
amplitudes = []
ratios     = []

export_activity = DataFrame();
data_groundtruth = npzread("data/ABidentities_blobel2021_chr2_35Mb_60Mb.npy")
mask = Vector{Bool}(data_groundtruth)

for i in 1:10
    name = names[i]
    file = h5open(["data/", name, "_discrete_n=3", ".hdf5"] |> join, "r")
    mat  = read(file["optimization/matrix"])
    vec  = read(file["optimization/vector"])
    offset = read(file["optimization/offset"])
    close(file)
    
    # define loss and gradient
    function loss(activity)
        transpose(activity) * mat * activity - 2*transpose(vec) * activity + offset
    end

    function loss_grad!(G, activity)
        G .= mat * activity .- vec
    end
    
    N = size(mat,1)
    start    = fill(0.001, N)
    lower    = fill(0.0, N)
    upper    = fill(Inf, N)

    result   = optimize(loss, loss_grad!, lower, upper, start, LBFGS() |> Fminbox,
            Optim.Options(show_trace=false)
    )
    
    activity = result.minimizer / mean(result.minimizer);
    
    append!(amplitudes, 
        2( mean(activity[mask]) - mean(activity[.!mask]) ) / ( mean(activity[mask]) + mean(activity[.!mask]) )
    )
    append!(ratios, 
        mean(activity[mask]) / mean(activity[.!mask])
    )
    
    export_activity[!, name] = activity;
    
    # plot profile
    #imshow(reshape(activity, 1, length(activity)), cmap=:coolwarm, extent=[1,1000,100,1], clim=[activitymin, activitymax]);
    imshow(reshape(activity, 1, length(activity)), cmap=:coolwarm, extent=[1,1000,100,1], clim=[0, 2]);
    xlim(left=1,right=1000);
    ylim(bottom=1,top=100);
    axis("off");
    tight_layout();
    savefig([folder, "/inferred_activity-", analysis_type, "-", name, ".png"] |> join,
        tight_layout=true, bbox_inches="tight", pad_inches=0, dpi=600);
    
    # print mean squared separation map of original data
    R, ΔR = ActivePolymer.Optimization.Interface.load_data(name);
    imshow(ΔR, cmap=cmasher.ocean, extent=[1,1000,1000,1], clim=[0,maximum(ΔR)]);
    xlim(left=1,right=1000);
    ylim(bottom=1,top=1000);
    axis("off");
    tight_layout();
    savefig([folder, "/inferred_activity-", analysis_type, "-original_separation-", name, ".png"] |> join,
        tight_layout=true, bbox_inches="tight", pad_inches=0, dpi=300);
    
    # print predicted mean squared separation map
    prediction = ActivePolymer.Optimization.Model.numeric_dense(result.minimizer, jacobian)
    imshow(prediction, cmap=cmasher.ocean, extent=[1,1000,1000,1], clim=[0,maximum(ΔR)]);
    xlim(left=1,right=1000);
    ylim(bottom=1,top=1000);
    axis("off");
    tight_layout();
    savefig([folder, "/inferred_activity-", analysis_type, "-predicted_separation-", name, ".png"] |> join,
        tight_layout=true, bbox_inches="tight", pad_inches=0, dpi=300);
end

data = DataFrame(ratio = values, delta = mappedvalues, deltaA = amplitudes, ratioA = ratios);
CSV.write([folder, "/inferred_activity-", analysis_type, "-amplitudes.csv"] |> join, data);

imshow(reshape(data_groundtruth, 1, length(data_groundtruth)), cmap=:coolwarm, extent=[1,1000,100,1], clim=[0, 1])
xlim(left=1,right=1000)
ylim(bottom=1,top=100)
axis("off")
tight_layout()
savefig([folder, "/ground_truth.png"] |> join,
    tight_layout=true, bbox_inches="tight", pad_inches=0, dpi=600)

export_activity[!, "position"] = Vector((1:size(export_activity,1)));
export_activity[!, "ground_truth"] = data_groundtruth;
CSV.write([folder, "/inferred_activity-", analysis_type, "-activity_vectors.csv"] |> join, export_activity);
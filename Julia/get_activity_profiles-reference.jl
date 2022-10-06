using NPZ
using Plots
using FFTW
using HDF5
using Optim
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

# include custom code
include("./Test/runtests.jl")
@everywhere include("./Source/ActivePolymer.jl")
using .ActivePolymer
using .ActivePolymer.CorrelationMatrices

# setup reference
reference = ActivePolymer.Optimization.DirectComparable.setup_reference_system(
    "Deq1", jacmodule=ActivePolymer.Jacobian.Discrete, modeltype=ActivePolymer.Optimization.Model.Full, n=3, padding=0.85)

# iterate over all parameters
for name in ["Deq1", "bcomps_2x", "bcomps_3x", "bcomps_4x", "bcomps_5x", "bcomps_7x", "bcomps_10x", "bcomps_19x", "bcomps_26x", "bcomps_39x"]
    ActivePolymer.Optimization.DirectComparable.setup_direct_system(name, reference, overwrite=true);
end
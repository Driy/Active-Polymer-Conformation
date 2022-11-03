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
using LinearAlgebra
using BenchmarkTools
using DelimitedFiles
using ImageFiltering

# include custom code
include("./Test/runtests.jl")
include("./Source/ActivePolymer.jl")
using .ActivePolymer
using .ActivePolymer.CorrelationMatrices

ΔRs = DataFrame()
Js  = DataFrame()
qs  = DataFrame()

for name in ["Deq1", "bcomps_2x", "bcomps_3x", "bcomps_4x", "bcomps_5x", "bcomps_7x", "bcomps_10x", "bcomps_19x", "bcomps_26x", "bcomps_39x"]
    R, ΔR = ActivePolymer.Optimization.Interface.load_data(name);
    ΔR_marginalized = ActivePolymer.Methods.Real.marginalize_translation(ΔR);
    N = size(ΔR_marginalized, 1)
    
    ΔRs[!, name] = ΔR_marginalized;
    Js[!, name]  = -1 ./ (sqrt(N) * dct(ΔR_marginalized)) |> x->imfilter(x, Kernel.gaussian((0,)));
    qs[!, name]  = [ActivePolymer.Methods.FastFourier._frequency_dct(i, N) for i in 1:N];    
end

# We fit the mechanical properties of a homogeneous polymer in thermal equilibrium
R, ΔR = ActivePolymer.Optimization.Interface.load_data("Deq1");

# Heuristically fix some properties of the polymer that cannot be fitted
jacmodule = ActivePolymer.Jacobian.Discrete;
n = 3;

# Polymer length
N = size(ΔR,1);

# Three-parameter fits
mech_pars = ActivePolymer.Optimization.Interface.fit_mechanics(
    ΔR, modeltype=ActivePolymer.Optimization.Model.Full, jacmodule=jacmodule, n=n, padding=0.85);

#
qfit          = [ActivePolymer.Methods.FastFourier._frequency_dct(i, N) for i in 1:N];    
qs[!, "fit"]  = qfit;    

#
model(q, p)   = (1/p[1])*jacmodule.J(p[2], p[3], n)(q);
Js[!, "fit"]  = model(qfit, mech_pars.minimizer);

#
ΔRs[!, "fit"] = ActivePolymer.Optimization.Model.numeric_marginalized(N, mech_pars.minimizer..., n=n, jacmodule=jacmodule);

# print fitted parameters
print(mech_pars.minimizer)

# export raw and fitted Jacobian 
mkpath("mechanics_fit")

for name in ["fit", "Deq1", "bcomps_2x", "bcomps_3x", "bcomps_4x", "bcomps_5x", "bcomps_7x", "bcomps_10x", "bcomps_19x", "bcomps_26x", "bcomps_39x"]
    data = DataFrame(q = qs[!, name], J = Js[!, name], deltaR = ΔRs[!, name], s = Vector(0:size(ΔRs[!, name],1)-1));
    CSV.write(["mechanics_fit/",name,".csv"] |> join, data);
end
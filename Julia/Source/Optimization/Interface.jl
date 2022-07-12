module Interface

using NPZ
using DelimitedFiles
using LinearAlgebra
using Statistics
using Optim

using ..Model
using ..Residual
using ...Methods

"""
Specific error for not invalid or not implemented transform types.
"""
struct InvalidType <:Exception
    msg::String
end

function load_data(name, path="../Share/June03/")
    name1 = ["comps_", name, "_rsquared"] |> join;
    name2 = ["mean_squared_separation_comp", name] |> join;
    R  = [path, name1, ".npy"] |> join |> NPZ.npzread;
    ΔR = [path, name2, ".csv"] |> join |> (x->readdlm(x, ',')) |> (x->x[begin+1:end, begin+1:end]) |> Matrix{Float64};
    return R, ΔR;
end

function fit_mechanics(ΔR; ModelType=Full, kwargs...)    
    # Three-parameter fits
    if ModelType==Model.Full
        parameters = [1.0, 0.001, 0.001]
        lower = [0.0, 0.0, 0.0]
        upper = [Inf, Inf, Inf]
    # Two-parameter fits
    elseif ModelType==Model.Reduced
        parameters = [1.0, 0.001]
        lower = [0.0, 0.0]
        upper = [Inf, Inf]
    else
        InvalidType("This Model Type is not defined!")
    end
    
    ΔR_marginalized = Methods.Real.marginalize_translation(ΔR);
    return optimize(
        Residual.numeric(ΔR_marginalized; kwargs...), 
        Residual.numeric_grad(ΔR_marginalized; kwargs...),
        lower, upper, parameters)
end

end
module Interface

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

function load_data(name, path="../Share/")
    R  = [path, "rsquared_", name, ".csv"] |> join |> (x->readdlm(x, ',', skipstart=1)) |> x->x[:] |> Vector{Float64};
    ΔR = [path, "mean_squared_separation_", name, ".csv"] |> join |> (x->readdlm(x, ',', skipstart=1)) |> Matrix{Float64};
    return R, ΔR;
end

function fit_mechanics(ΔR; modeltype=Model.Full, kwargs...)    
    # Three-parameter fits
    if modeltype==Model.Full
        parameters = [1.0, 0.001, 0.001]
        lower = [0.0, 0.0, 0.0]
        upper = [Inf, Inf, Inf]
    # Two-parameter fits
    elseif modeltype==Model.Reduced
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

function extract_mechanics(ΔR; window=1)
    ΔR_marginalized = ActivePolymer.Methods.Real.marginalize_translation(ΔR);
    N = size(ΔR_marginalized,1);
    J = -1 ./ (sqrt(N) * dct(ΔR_marginalized)) |> x->imfilter(x, Kernel.gaussian((window,)));
    return J;
end

end
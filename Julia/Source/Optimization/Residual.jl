module Residual

using LinearAlgebra
using Statistics

using ..Model
using ...Jacobian
using ...CorrelationMatrices

"""
numeric(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Determine the mean squared error of a proposed model when compared to the marginalized mean squared separation data `ΔR_marginalized`.
"""
function numeric(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function (parameters)
        residual_vector = Model.numeric_marginalized(N, parameters...; kwargs...) - ΔR_marginalized;
        squared_error_vector = residual_vector.^2;
        return squared_error_vector[begin:end-N_padding] |> mean;
    end
end

"""
numeric(ΔR::AbstractMatrix, parameters; padding::Real=0.0, kwargs...)

Determine the mean squared error of a proposed model when compared to the mean squared separation data `ΔR_marginalized`.
"""
function numeric(ΔR::AbstractMatrix, jacobian; padding::Real=0.0)
    N, N_padding = (1, padding) .* size(ΔR, 1) .|> Int64;
    return function(activity::AbstractVector)
        residual_matrix = Model.numeric_dense(activity, jacobian) - ΔR;
        squared_error_matrix = residual_matrix.^2;
        return squared_error_matrix[begin+N_padding:end-N_padding, begin+N_padding:end-N_padding] |> mean;
    end
end

"""
numeric_grad(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Method that determines the gradient with respect to the `parameters`, of the mean squared error of a proposed model when compared to the marginalized mean squared separation data `ΔR_marginalized`.
"""
function numeric_grad(ΔR_marginalized::AbstractVector; padding::Real=0.75, jacmodule=Jacobian.Standard, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function(G, parameters)
        residual_vector = Model.numeric_marginalized(N, parameters...; jacmodule = jacmodule, kwargs...) - ΔR_marginalized;
        grad_T_val = 2residual_vector .* Model.numeric_marginalized_grad(N, parameters...; jacmodule = jacmodule, kwargs...);
        grad_α_val = 2residual_vector .* Model.numeric_marginalized_grad(jacmodule.dJ_dα, N, parameters...; jacmodule = jacmodule, kwargs...);
        grad_κ_val = 2residual_vector .* Model.numeric_marginalized_grad(jacmodule.dJ_dκ, N, parameters...; jacmodule = jacmodule, kwargs...);
        result_vector = [grad_vector[begin:end-N_padding] |> mean for grad_vector in [grad_T_val, grad_α_val, grad_κ_val]];
        
        G[:] .= result_vector[begin:size(parameters,1)];
    end
end

"""
numeric_grad(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Method that determines the gradient with respect to the local activity, of the mean squared error of a proposed model when compared to the mean squared separation data `ΔR`.
"""
function numeric_grad(ΔR::AbstractMatrix, jacobian; padding::Real=0.0)
    N, N_padding = (1, padding) .* size(ΔR, 1) .|> Int64;
    return function(G, activity)
        residual_matrix = Model.numeric_dense(activity, jacobian) - ΔR;
        for i in 1:N
            G[i] = 2residual_matrix .* Model.numeric_dense(
                CorrelationMatrices.diagonal_delta(i, N), jacobian) |> mean;
        end
    end    
end

"""
numeric_grad_bruteforce(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Brute force method.

Method that determines the gradient with respect to the local activity, of the mean squared error of a proposed model when compared to the mean squared separation data `ΔR`.
"""
function numeric_grad_bruteforce(ΔR::AbstractMatrix, jacobian; padding::Real=0.0, tolerance=1e-5)
    N, N_padding = (1, padding) .* size(ΔR, 1) .|> Int64;
    tmpfun = numeric(ΔR, jacobian; padding=padding)
    return function(G, activity)
        for i in 1:N
            d_activity = zeros(N); d_activity[i] = tolerance;
            G[i] = (tmpfun(activity + d_activity) - tmpfun(activity - d_activity)) / 2tolerance;            
        end
    end    
end

"""
analytic(ΔR_marginalized::AbstractVector; padding::Real=0.75)

Determine the mean squared error of a proposed analytic model when compared to the marginalized mean squared separation data `ΔR_marginalized`.
"""
function analytic(ΔR_marginalized::AbstractVector; padding::Real=0.75)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function(parameters)
        residual_vector = Model.analytic_marginalized(N, parameters...) - ΔR_marginalized;
        squared_error_vector = residual_vector.^2;
        return squared_error_vector[begin:end-N_padding] |> mean;
    end
end

end
module MethodsFitting

using LinearAlgebra
using Statistics

using ..TransformForward
using ..StandardFunctions
using ..MethodsAnalytic
using ..MethodsReal
using ..WrapperFFTW

function model_numeric_dense(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)
    return TransformForward.compute_conformation(
        T, N, J=(q->StandardFunctions.J(q; α=α, κ=κ, n=n)), fourier_type=WrapperFFTW.DCT) |> MethodsReal.correlation_to_separation;
end

function model_numeric_dense(T::AbstractVector, α::Real=0, κ::Real=0; n::Real=4)
    # compute semianalytic solution for homogeneous activity matrix
    return TransformForward.compute_conformation(
        T, J=(q->StandardFunctions.J(q; α=α, κ=κ, n=n)), fourier_type=WrapperFFTW.DCT) |> MethodsReal.correlation_to_separation;
end

function model_numeric_dense_grad(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)
    return TransformForward.compute_conformation(
        1.0, N, J=(q->StandardFunctions.J(q; α=α, κ=κ, n=n)), fourier_type=WrapperFFTW.DCT) |> MethodsReal.correlation_to_separation;
end

function model_numeric_dense_grad(dJ_dα::Function, N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)
    return TransformForward.compute_conformation_grad(
        T, N, J=(q->StandardFunctions.J(q; α=α, κ=κ, n=n)), dJ_dα=(q->dJ_dα(q; α=α, κ=κ, n=n)), 
        fourier_type=WrapperFFTW.DCT) |> MethodsReal.correlation_to_separation;
end

function model_numeric_marginalized(parameters...; kwargs...)
    return model_numeric_dense(parameters...; kwargs...) |> MethodsReal.marginalize_translation;
end

function model_numeric_marginalized_grad(parameters...; kwargs...)
    return model_numeric_dense_grad(parameters...; kwargs...) |> MethodsReal.marginalize_translation;
end

function model_numeric_marginalized_grad(dJ_dα::Function, parameters...; kwargs...)
    return model_numeric_dense_grad(dJ_dα, parameters...; kwargs...) |> MethodsReal.marginalize_translation;
end

function residual_numeric(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function (parameters)
        residual_vector = model_numeric_marginalized(N, parameters...; kwargs...) - ΔR_marginalized;
        squared_error_vector = residual_vector.^2;
        return squared_error_vector[begin:end-N_padding] |> mean;
    end
end

function residual_numeric(ΔR::AbstractMatrix, parameters; padding::Real=0.0, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR, 1) .|> Int64;
    return function(activity::AbstractVector)
        residual_matrix = model_numeric_dense(activity, parameters[2:end]...; kwargs...) - ΔR;
        squared_error_matrix = residual_matrix.^2;
        return squared_error_matrix[begin+N_padding:end-N_padding, begin+N_padding:end-N_padding] |> mean;
    end
end

function residual_numeric_grad(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function(G, parameters)
        residual_vector = model_numeric_marginalized(N, parameters...; kwargs...) - ΔR_marginalized;
        grad_T_val = 2residual_vector .* model_numeric_marginalized_grad(N, parameters...);
        grad_α_val = 2residual_vector .* model_numeric_marginalized_grad(StandardFunctions.dJ_dα, N, parameters...);
        grad_κ_val = 2residual_vector .* model_numeric_marginalized_grad(StandardFunctions.dJ_dκ, N, parameters...);
        result_vector = [grad_vector[begin:end-N_padding] |> mean for grad_vector in [grad_T_val, grad_α_val, grad_κ_val]];
        
        G[:] .= result_vector[begin:size(parameters,1)];        
    end
end

function model_analytic_marginalized(Δs::AbstractVector, T::Real, α::Real, κ::Real)
    return MethodsAnalytic.separation_generic(Δs; T=T, α=α, κ=κ);    
end

function residual_analytic(ΔR_marginalized::AbstractVector; padding::Real=0.75)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    Δs = Vector{Float64}(0:N-1);
    return function(parameters)
        residual_vector = model_analytic_marginalized(Δs, parameters...) - ΔR_marginalized;
        squared_error_vector = residual_vector.^2;
        return squared_error_vector[begin:end-N_padding] |> mean;
    end
end

end
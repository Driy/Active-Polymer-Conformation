module MethodsFitting

using LinearAlgebra
using Statistics

using ..CorrelationMatrices

using ..Transform
using ..Jacobian
using ..Methods
using ..Methods.FastFourier

"""
model_numeric_dense(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)

Compute mean squared separation between different monomers of an `N`-Polymer with a confining harmonic spring `α`, bending rigidity `κ`, and a homogeneous level of activity `T`.
"""
function model_numeric_dense(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    return Transform.Forward.compute_conformation(
        T, N, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
model_numeric_dense(T::AbstractVector, α::Real=0, κ::Real=0; n::Real=4)

Compute mean squared separation between different monomers of a Polymer with a confining harmonic spring `α`, bending rigidity `κ`, and an inhomogeneous level of activity given by the vector `T`.
"""
function model_numeric_dense(T::AbstractVector, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    # compute semianalytic solution for homogeneous activity matrix
    return Transform.Forward.compute_conformation(
        T, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
model_numeric_dense(T::AbstractMatrix, α::Real=0, κ::Real=0; n::Real=4)

Compute mean squared separation between different monomers of a Polymer with a confining harmonic spring `α`, bending rigidity `κ`, and the activity given by the matrix `T`.
"""
function model_numeric_dense(T::AbstractMatrix, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    # compute semianalytic solution for homogeneous activity matrix
    return Transform.Forward.compute_conformation(
        T, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
model_numeric_dense_grad(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)

Determine how a homogeneous change in the activity affects the mean squared separation between different monomers.
"""
function model_numeric_dense_grad(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    return Transform.Forward.compute_conformation(
        1.0, N, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
model_numeric_dense_grad(dJ_dα::Function, N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)

Determine how a change in the mechanical properties of the polymer, encoded by the change in Jacobian `dJ_dα`, affects the mean squared separation between different monomers.
"""
function model_numeric_dense_grad(dJ_dα::Function, N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    return Transform.Forward.compute_conformation_grad(
        T, N, J=jacmodule.J(α, κ, n), dJ_dα=dJ_dα(α, κ, n), 
        fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
model_numeric_marginalized(parameters...; kwargs...)

Marginalized mean squared separation between different monomers of a polymer, for a translationally invariant system where the mean squared separation only depends on the distance in sequence space.
"""
function model_numeric_marginalized(parameters...; kwargs...)
    return model_numeric_dense(parameters...; kwargs...) |> Methods.Real.marginalize_translation;
end

"""
model_numeric_marginalized_grad(parameters...; kwargs...)

Determine how a homogeneous change in the activity affects the (marginalized) mean squared separation between different monomers.
"""
function model_numeric_marginalized_grad(parameters...; kwargs...)
    return model_numeric_dense_grad(parameters...; kwargs...) |> Methods.Real.marginalize_translation;
end

"""
model_numeric_marginalized_grad(dJ_dα::Function, parameters...; kwargs...)

Determine how a change in the mechanical properties of the polymer, encoded by the change in Jacobian `dJ_dα`, affects the (marginalized) mean squared separation between different monomers.
"""
function model_numeric_marginalized_grad(dJ_dα::Function, parameters...; kwargs...)
    return model_numeric_dense_grad(dJ_dα, parameters...; kwargs...) |> Methods.Real.marginalize_translation;
end

"""
residual_numeric(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Determine the mean squared error of a proposed model when compared to the marginalized mean squared separation data `ΔR_marginalized`.
"""
function residual_numeric(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function (parameters)
        residual_vector = model_numeric_marginalized(N, parameters...; kwargs...) - ΔR_marginalized;
        squared_error_vector = residual_vector.^2;
        return squared_error_vector[begin:end-N_padding] |> mean;
    end
end

"""
residual_numeric(ΔR::AbstractMatrix, parameters; padding::Real=0.0, kwargs...)

Determine the mean squared error of a proposed model when compared to the mean squared separation data `ΔR_marginalized`.
"""
function residual_numeric(ΔR::AbstractMatrix, parameters; padding::Real=0.0, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR, 1) .|> Int64;
    return function(activity::AbstractVector)
        residual_matrix = model_numeric_dense(activity, parameters[2:end]...; kwargs...) - ΔR;
        squared_error_matrix = residual_matrix.^2;
        return squared_error_matrix[begin+N_padding:end-N_padding, begin+N_padding:end-N_padding] |> mean;
    end
end

"""
residual_numeric_grad(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Method that determines the gradient with respect to the `parameters`, of the mean squared error of a proposed model when compared to the marginalized mean squared separation data `ΔR_marginalized`.
"""
function residual_numeric_grad(ΔR_marginalized::AbstractVector; padding::Real=0.75, jacmodule=Jacobian.Standard, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function(G, parameters)
        residual_vector = model_numeric_marginalized(N, parameters...; jacmodule = jacmodule, kwargs...) - ΔR_marginalized;
        grad_T_val = 2residual_vector .* model_numeric_marginalized_grad(N, parameters...; jacmodule = jacmodule, kwargs...);
        grad_α_val = 2residual_vector .* model_numeric_marginalized_grad(jacmodule.dJ_dα, N, parameters...; jacmodule = jacmodule, kwargs...);
        grad_κ_val = 2residual_vector .* model_numeric_marginalized_grad(jacmodule.dJ_dκ, N, parameters...; jacmodule = jacmodule, kwargs...);
        result_vector = [grad_vector[begin:end-N_padding] |> mean for grad_vector in [grad_T_val, grad_α_val, grad_κ_val]];
        
        G[:] .= result_vector[begin:size(parameters,1)];
    end
end

"""
residual_numeric_grad(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Method that determines the gradient with respect to the local activity, of the mean squared error of a proposed model when compared to the mean squared separation data `ΔR`.
"""
function residual_numeric_grad(ΔR::AbstractMatrix, parameters; padding::Real=0.0, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR, 1) .|> Int64;
    return function(G, activity)
        residual_matrix = model_numeric_dense(activity, parameters[2:end]...; kwargs...) - ΔR;
        for i in 1:N
            G[i] = 2residual_matrix .* model_numeric_dense(
                CorrelationMatrices.diagonal_delta(i, N), parameters[2:end]...; kwargs...) |> mean;
        end
    end    
end

"""
residual_numeric_grad_bruteforce(ΔR_marginalized::AbstractVector; padding::Real=0.75, kwargs...)

Brute force method.

Method that determines the gradient with respect to the local activity, of the mean squared error of a proposed model when compared to the mean squared separation data `ΔR`.
"""
function residual_numeric_grad_bruteforce(ΔR::AbstractMatrix, parameters; padding::Real=0.0, tolerance=1e-5, kwargs...)
    N, N_padding = (1, padding) .* size(ΔR, 1) .|> Int64;
    tmpfun = residual_numeric(ΔR, parameters; padding=padding, kwargs...)
    return function(G, activity)
        for i in 1:N
            d_activity = zeros(N); d_activity[i] = tolerance;
            G[i] = (tmpfun(activity + d_activity) - tmpfun(activity - d_activity)) / 2tolerance;            
        end
    end    
end

"""
model_analytic_marginalized(Δs::AbstractVector, T::Real, α::Real, κ::Real)

Analytic expression for the marginalized mean squared separation between different monomers of a continuous polymer, for a translationally invariant system where the mean squared separation only depends on the distance in sequence space.
"""
function model_analytic_marginalized(N::Integer, T::Real, α::Real, κ::Real)
    Δs = Vector{Float64}(0:N-1);
    return Methods.Analytic.separation_generic(Δs; T=T, α=α, κ=κ);    
end

"""
residual_analytic(ΔR_marginalized::AbstractVector; padding::Real=0.75)

Determine the mean squared error of a proposed analytic model when compared to the marginalized mean squared separation data `ΔR_marginalized`.
"""
function residual_analytic(ΔR_marginalized::AbstractVector; padding::Real=0.75)
    N, N_padding = (1, padding) .* size(ΔR_marginalized, 1) .|> Int64;
    return function(parameters)
        residual_vector = model_analytic_marginalized(N, parameters...) - ΔR_marginalized;
        squared_error_vector = residual_vector.^2;
        return squared_error_vector[begin:end-N_padding] |> mean;
    end
end

"""
coupling_matrix!(mat, i, j, parameters; kwargs...)

Determine elements of the coupling matrix, which determines the gradient of the [proposed mean squared separation matrix] squared, with respect to localized changes in activity.
"""
function populate_coupling_matrix!(mat, i, j, parameters; kwargs...)
    N = size(mat, 1)
    tmp1 = model_numeric_dense(
        CorrelationMatrices.diagonal_delta(i, N), parameters[2:end]...; kwargs...);
    tmp2 = model_numeric_dense(
        CorrelationMatrices.diagonal_delta(j, N), parameters[2:end]...; kwargs...);
    mat[i,j] = mean(tmp1 .* tmp2);
end

"""
coupling_vector!(mat, i, j, parameters; kwargs...)

Determine elements of the coupling vector, which determines the gradient of the [proposed mean squared separation matrix] times the [mean squared separation data], with respect to localized changes in activity.
"""
function populate_coupling_vector!(vec, i, ΔR, parameters; kwargs...)
    N = size(vec, 1)
    tmp1 = model_numeric_dense(
        CorrelationMatrices.diagonal_delta(i, N), parameters[2:end]...; kwargs...);
    vec[i] = mean(tmp1 .* ΔR);
end


end
module Model

using LinearAlgebra
using Statistics

using ...CorrelationMatrices

using ...Transform
using ...Jacobian
using ...Methods
using ...Methods.FastFourier

@enum Type Full=0 Reduced=1

"""
numeric_dense(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)

Compute mean squared separation between different monomers of an `N`-Polymer with a confining harmonic spring `α`, bending rigidity `κ`, and a homogeneous level of activity `T`.
"""
function numeric_dense(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    return Transform.Forward.compute_conformation(
        T, N, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
numeric_dense(T::AbstractVector, α::Real=0, κ::Real=0; n::Real=4)

Compute mean squared separation between different monomers of a Polymer with a confining harmonic spring `α`, bending rigidity `κ`, and an inhomogeneous level of activity given by the vector `T`.
"""
function numeric_dense(T::AbstractVector, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    # compute semianalytic solution for homogeneous activity matrix
    return Transform.Forward.compute_conformation(
        T, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
numeric_dense(T::AbstractMatrix, α::Real=0, κ::Real=0; n::Real=4)

Compute mean squared separation between different monomers of a Polymer with a confining harmonic spring `α`, bending rigidity `κ`, and the activity given by the matrix `T`.
"""
function numeric_dense(T::AbstractMatrix, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    # compute semianalytic solution for homogeneous activity matrix
    return Transform.Forward.compute_conformation(
        T, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
numeric_dense_grad(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)

Determine how a homogeneous change in the activity affects the mean squared separation between different monomers.
"""
function numeric_dense_grad(N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    return Transform.Forward.compute_conformation(
        1.0, N, J=jacmodule.J(α, κ, n), fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
numeric_dense_grad(dJ_dα::Function, N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4)

Determine how a change in the mechanical properties of the polymer, encoded by the change in Jacobian `dJ_dα`, affects the mean squared separation between different monomers.
"""
function numeric_dense_grad(dJ_dα::Function, N::Integer, T::Real, α::Real=0, κ::Real=0; n::Real=4, jacmodule=Jacobian.Standard)
    return Transform.Forward.compute_conformation_grad(
        T, N, J=jacmodule.J(α, κ, n), dJ_dα=dJ_dα(α, κ, n), 
        fourier_type=FastFourier.DCT) |> Methods.Real.correlation_to_separation;
end

"""
numeric_marginalized(parameters...; kwargs...)

Marginalized mean squared separation between different monomers of a polymer, for a translationally invariant system where the mean squared separation only depends on the distance in sequence space.
"""
function numeric_marginalized(parameters...; kwargs...)
    return numeric_dense(parameters...; kwargs...) |> Methods.Real.marginalize_translation;
end

"""
numeric_marginalized_grad(parameters...; kwargs...)

Determine how a homogeneous change in the activity affects the (marginalized) mean squared separation between different monomers.
"""
function numeric_marginalized_grad(parameters...; kwargs...)
    return numeric_dense_grad(parameters...; kwargs...) |> Methods.Real.marginalize_translation;
end

"""
numeric_marginalized_grad(dJ_dα::Function, parameters...; kwargs...)

Determine how a change in the mechanical properties of the polymer, encoded by the change in Jacobian `dJ_dα`, affects the (marginalized) mean squared separation between different monomers.
"""
function numeric_marginalized_grad(dJ_dα::Function, parameters...; kwargs...)
    return numeric_dense_grad(dJ_dα, parameters...; kwargs...) |> Methods.Real.marginalize_translation;
end

"""
analytic_marginalized(Δs::AbstractVector, T::Real, α::Real, κ::Real)

Analytic expression for the marginalized mean squared separation between different monomers of a continuous polymer, for a translationally invariant system where the mean squared separation only depends on the distance in sequence space.
"""
function analytic_marginalized(N::Integer, T::Real, α::Real, κ::Real)
    Δs = Vector{Float64}(0:N-1);
    return Methods.Analytic.separation_generic(Δs; T=T, α=α, κ=κ);    
end

end
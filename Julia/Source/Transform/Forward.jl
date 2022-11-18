module Forward

using LinearAlgebra

using ...Jacobian
using ...Methods
using ...Methods.FastFourier

export compute_conformation, compute_conformation_grad

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(matrix_activity::AbstractMatrix, J::Function = Jacobian.Standard.J₀;
        fourier_type = FastFourier.DCT)
    # transform to Fourier space; always dense!
    tmp = FastFourier.forward(matrix_activity, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    Methods.Spectral.activity_to_correlation!(tmp, J, fourier_type = fourier_type);
    # transform to real space
    FastFourier.backward!(tmp, fourier_type = fourier_type);
    # return mean square separation
    return tmp;
end

"""
compute_mean_squared_traveled_distance(matrix_activity::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Calculate the mean squared traveled distance after the time τ, for all loci along a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_mean_squared_traveled_distance(matrix_activity::AbstractMatrix, J::Function, τ::Real;
        fourier_type = FastFourier.DCT)
    # transform to Fourier space; always dense!
    tmp = FastFourier.forward(matrix_activity, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    Methods.Spectral.activity_to_correlation!(tmp, J, τ, fourier_type = fourier_type);
    # transform to real space
    FastFourier.backward!(tmp, fourier_type = fourier_type);
    # return mean square separation
    return 2tmp |> diag;
end

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(matrix_activity::AbstractMatrix, J::AbstractVector;
        fourier_type = FastFourier.DCT)
    # transform to Fourier space; always dense!
    tmp = FastFourier.forward(matrix_activity, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    Methods.Spectral.activity_to_correlation!(tmp, J, fourier_type = fourier_type);
    # transform to real space
    FastFourier.backward!(tmp, fourier_type = fourier_type);
    # return mean square separation
    return tmp;
end

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function, dJ_dα::Function, fourier_type = FastFourier.DCT)

This is a gradient method that computes the response to changing a parameter in the Jacobian. Takes the derivative of the Jacobian as keyword argument.

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(matrix_activity::AbstractMatrix, J::Function, dJ_dα::Function; 
        fourier_type = FastFourier.DCT)
    # transform to Fourier space; always dense!
    tmp = FastFourier.forward(matrix_activity, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    Methods.Spectral.activity_to_correlation!(tmp, J, dJ_dα, fourier_type = fourier_type);
    # transform to real space
    FastFourier.backward!(tmp, fourier_type = fourier_type);
    # return mean square separation
    return tmp;
end

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(vector_activity::AbstractVector, args...; kwargs...)
    return compute_conformation(vector_activity |> diagm, args...; kwargs...);
end

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(scalar_activity::Real, N::Integer, args...; kwargs...)
    return compute_conformation(fill(scalar_activity, N) |> diagm, args...; kwargs...);
end

end
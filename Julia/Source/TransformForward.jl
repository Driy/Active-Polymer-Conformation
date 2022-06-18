module TransformForward

using LinearAlgebra

using ..WrapperFFTW
using ..MethodsReal
using ..MethodsSpectral
using ..StandardFunctions

export compute_conformation

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(matrix_activity::AbstractMatrix; 
        J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)
    # transform to Fourier space; always dense!
    tmp = WrapperFFTW.forward(matrix_activity, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    MethodsSpectral.activity_to_correlation!(tmp, J = J, fourier_type = fourier_type);
    # transform to real space
    WrapperFFTW.backward!(tmp, fourier_type = fourier_type);
    # return mean square separation
    return tmp;
end

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(vector_activity::AbstractVector; 
        J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)
    return compute_conformation(vector_activity |> diagm, J=J, fourier_type=fourier_type);
end

"""
compute_conformation(matrix_activity::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)

Calculate the mean square separation of a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_conformation(scalar_activity::Float64, N::Int64; 
        J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)
    return compute_conformation(fill(scalar_activity, N) |> diagm, J=J, fourier_type=fourier_type);
end

end
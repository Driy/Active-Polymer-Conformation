module TransformForward

using ..WrapperFFTW
using ..MethodsReal
using ..MethodsSpectral
using ..StandardFunctions

export compute_conformation

"""
mean_square_separation(matrix::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)

Calculate the mean square separation of a polymer with given correlation matrix C and Jacobian J.
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
    return MethodsReal.correlation_to_position(tmp), MethodsReal.correlation_to_separation(tmp);
end

end
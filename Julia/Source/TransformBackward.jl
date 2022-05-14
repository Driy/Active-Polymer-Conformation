module TransformBackward

using ..WrapperFFTW
using ..MethodsReal
using ..MethodsSpectral
using ..StandardFunctions

export extract_excitations

"""
extract_excitations(vector_position::AbstractVector, matrix_separation::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)

Extract correlation function from mean square separation `matrix`, for a polymer with Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. 
"""
function extract_excitations(vector_position::AbstractVector, matrix_separation::AbstractMatrix; 
        J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)
    
    tmp = MethodsReal.position_and_separation_to_correlation(vector_position, matrix_separation)
    # transform to Fourier space; always dense!
    tmp = WrapperFFTW.forward(tmp, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    MethodsSpectral.correlation_to_activity!(tmp, J = J, fourier_type = fourier_type);
    # transform to real space
    WrapperFFTW.backward!(tmp, fourier_type = fourier_type);
    return tmp |> real;
end

"""
extract_excitations(matrix_separation::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)

Extract correlation function from mean square separation `matrix`, for a polymer with Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. 
"""
function extract_excitations(matrix_separation::AbstractMatrix; 
        J::Function = StandardFunctions.J₀, fourier_type = WrapperFFTW.DCT)
    
    tmp = MethodsReal.position_and_separation_to_correlation(
        zeros(size(matrix_separation)[1]), 
        matrix_separation)
    # transform to Fourier space; always dense!
    tmp = WrapperFFTW.forward(tmp, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    MethodsSpectral.correlation_to_activity!(tmp, J = J, fourier_type = fourier_type);
    # transform to real space
    WrapperFFTW.backward!(tmp, fourier_type = fourier_type);
    return tmp |> real;
end

end
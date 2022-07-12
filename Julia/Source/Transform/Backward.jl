module Backward

using ...Jacobian
using ...Methods
using ...Methods.FastFourier

export extract_excitations

"""
extract_excitations(vector_position::AbstractVector, matrix_separation::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Extract correlation function from mean square separation `matrix`, for a polymer with Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. 
"""
function extract_excitations(vector_position::AbstractVector, matrix_separation::AbstractMatrix; 
        J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)
    
    tmp = Methods.Real.position_and_separation_to_correlation(vector_position, matrix_separation)
    # transform to Fourier space; always dense!
    tmp = FastFourier.forward(tmp, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    Methods.Spectral.correlation_to_activity!(tmp, J = J, fourier_type = fourier_type);
    # transform to real space
    FastFourier.backward!(tmp, fourier_type = fourier_type);
    return tmp |> real;
end

"""
extract_excitations(matrix_separation::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Extract correlation function from mean square separation `matrix`, for a polymer with Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. 
"""
function extract_excitations(matrix_separation::AbstractMatrix; 
        J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)
    
    tmp = Methods.Real.position_and_separation_to_correlation(
        zeros(size(matrix_separation)[1]), 
        matrix_separation)
    # transform to Fourier space; always dense!
    tmp = FastFourier.forward(tmp, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    Methods.Spectral.correlation_to_activity!(tmp, J = J, fourier_type = fourier_type);
    # transform to real space
    FastFourier.backward!(tmp, fourier_type = fourier_type);
    return tmp |> real;
end

end
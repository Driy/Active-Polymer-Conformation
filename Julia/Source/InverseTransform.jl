module InverseTransform

using ..FourierTransform
using ..StandardFunctions

export extract_excitations

"""
fourier_ΔR_to_C!(matrix::AbstractMatrix; J::Function, fourier_type)

Map the mean square separation in Fourier space `matrix` to return the corresponding noise correlation matrix C of the polymer. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. Specify `J` to define the Jacobian of the polymer. 
"""
function fourier_ΔR_to_C!(matrix::AbstractMatrix; J::Function, fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = FourierTransform.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            matrix[id] *= -0.5;
            matrix[id] *= J(q) + J(k);
        end
    end
    # the homogeneous modes cannot be determined, because ΔR is invariant with respect to them!
    matrix[begin,:] .= 0.0; matrix[:,begin] .= 0.0;
end

"""
extract_excitations(matrix::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = FourierTransform.DCT)

Extract correlation function from mean square separation `matrix`, for a polymer with Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. 
"""
function extract_excitations(matrix::AbstractMatrix; 
        J::Function = StandardFunctions.J₀, fourier_type = FourierTransform.DCT)
    # transform to Fourier space; always dense!
    tmp = FourierTransform.forward(matrix, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    fourier_ΔR_to_C!(tmp, J = J, fourier_type = fourier_type);
    # transform to real space
    FourierTransform.inverse!(tmp, fourier_type = fourier_type);
    return real(tmp);
end

end
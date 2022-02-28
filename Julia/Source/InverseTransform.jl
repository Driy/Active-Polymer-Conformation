module InverseTransform

using ..FourierTransform
using ..StandardFunctions

export extract_excitations

"""
fourier_ΔR_to_C!(fourier_type, matrix::AbstractMatrix, J::Function)

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
extract_excitations(fourier_type, matrix::AbstractMatrix, J::Function = J₀)

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

"""
test_extraction(C)

Test the extraction of the correlation function for a given correlation function.
Note that we cannot extract information about homogeneous contributions!
Thus, this ONLY works perfectly if there are no homogeneous contributions.
"""
function test_extraction(matrix::AbstractMatrix)
    R = mean_square_separation(matrix);
    C_extracted = extract_correlations(R);
    return (mean(matrix - C_extracted), std(matrix - C_extracted));
end

end
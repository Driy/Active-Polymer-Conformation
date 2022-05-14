module MethodsSpectral

using ..WrapperFFTW

export activity_to_correlation!, correlation_to_activity!

"""
activity_to_correlation!(matrix::AbstractMatrix; J::Function, fourier_type)

Map the noise correlation matrix `tmp` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

Since we are dealing with discrete Fourier transforms, it is tedious to extend this method to determine the tangent-tangent autocorrelation function. For that specific purpose, there is the separate discrete method [`real_R_to_ttacf`](@doc), which operates in real space.

## Example
```julia-repl
julia> fourier_C_to_R!(DCT, C, J₀);
```
"""
function activity_to_correlation!(matrix::AbstractMatrix; J::Function, fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = WrapperFFTW.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            matrix[id] /= J(q) + J(k);
        end
    end
    # Make sure that the homogeneous mode is not NaN! 
    # It is irrelevant for the mean square separation anyways.
    if !isfinite(matrix[begin,begin])
        matrix[begin,begin] = 0
    end
end

"""
correlation_to_activity!(matrix::AbstractMatrix; J::Function, fourier_type)

Map the position correlation in Fourier space `matrix` to return the corresponding noise correlation matrix C of the polymer. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. Specify `J` to define the Jacobian of the polymer. 
"""
function correlation_to_activity!(matrix::AbstractMatrix; J::Function, fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = WrapperFFTW.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            matrix[id] *= J(q) + J(k);
        end
    end
end

end
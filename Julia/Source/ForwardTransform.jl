module ForwardTransform

using ..FourierTransform
using ..StandardFunctions

export mean_square_separation

"""
fourier_C_to_R!(matrix::AbstractMatrix; J::Function, fourier_type)

Map the noise correlation matrix `tmp` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

Since we are dealing with discrete Fourier transforms, it is tedious to extend this method to determine the tangent-tangent autocorrelation function. For that specific purpose, there is the separate discrete method [`real_R_to_ttacf`](@doc), which operates in real space.

## Example
```julia-repl
julia> fourier_C_to_R!(DCT, C, J₀);
```
"""
function fourier_C_to_R!(matrix::AbstractMatrix; J::Function, fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = FourierTransform.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            matrix[id] /= J(q) + J(k);
        end
    end
    # Make sure that the homogeneous mode is not NaN! 
    # It is irrelevant for the mean square separation anyways.
    matrix[begin,:] .= 0.0; matrix[:,begin] .= 0.0;
end

"""
real_R_to_ΔR(matrix::AbstractMatrix)

Calculate mean square separation matrix from correlation matrix between different Rouse modes (`matrix`).

## Example
```julia-repl
julia> ΔR = real_R_to_ΔR(R);
```
"""
function real_R_to_ΔR(matrix::AbstractMatrix)
    return [
        let (x,y) = id.I
            matrix[x,x] + matrix[y,y] - 2matrix[x,y]
        end
        for id in CartesianIndices(matrix)];
end

"""
mean_square_separation(matrix::AbstractMatrix; J::Function = StandardFunctions.J₀, fourier_type = FourierTransform.DCT)

Calculate the mean square separation of a polymer with given correlation matrix C and Jacobian J.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function mean_square_separation(matrix::AbstractMatrix; 
        J::Function = StandardFunctions.J₀, fourier_type = FourierTransform.DCT)
    # transform to Fourier space; always dense!
    tmp = FourierTransform.forward(matrix, fourier_type = fourier_type);
    # manipulate as needed in Fourier space
    fourier_C_to_R!(tmp, J = J, fourier_type = fourier_type);
    # transform to real space
    FourierTransform.inverse!(tmp, fourier_type = fourier_type);
    # return mean square separation
    return real(real_R_to_ΔR(tmp));
end

end
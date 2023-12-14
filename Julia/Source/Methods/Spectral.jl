module Spectral

using ..FastFourier

export activity_to_correlation!, correlation_to_activity!

"""
activity_to_correlation!(matrix::AbstractMatrix; J::Function, fourier_type)

Map the noise correlation matrix `matrix` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with diagonal Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

Since we are dealing with discrete Fourier transforms, it is tedious to extend this method to determine the tangent-tangent autocorrelation function. For that specific purpose, there is the separate discrete method [`real_R_to_ttacf`](@doc), which operates in real space.

## Example
```julia-repl
julia> activity_to_correlation!(C, J=J₀, fourier_type=FastFourier.DCT);
```
"""
function activity_to_correlation!(matrix::AbstractMatrix, J::AbstractVector; fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = id.I
            # manipulate as needed in Fourier space
            local val = 1.0/(J[q] + J[k])
            # Make sure that the homogeneous mode is not NaN! 
            # It is irrelevant for the mean square separation anyways.
            if !isfinite(val) val = 0 end
            matrix[id] *= val;
        end
    end
end

"""
activity_to_correlation!(matrix::AbstractMatrix; J::Function, fourier_type)

Map the noise correlation matrix `matrix` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with diagonal Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

Since we are dealing with discrete Fourier transforms, it is tedious to extend this method to determine the tangent-tangent autocorrelation function. For that specific purpose, there is the separate discrete method [`real_R_to_ttacf`](@doc), which operates in real space.

## Example
```julia-repl
julia> activity_to_correlation!(C, J=J₀, fourier_type=FastFourier.DCT);
```
"""
function activity_to_correlation!(matrix::AbstractMatrix, J::Function; fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            local val = 1.0/(J(q) + J(k))
            # Make sure that the homogeneous mode is not NaN! 
            # It is irrelevant for the mean square separation anyways.
            if !isfinite(val) val = 0 end
            matrix[id] *= val;
        end
    end
end

"""
activity_to_correlation!(matrix::AbstractMatrix; J::Function, τ::Real, fourier_type)

Map the noise correlation matrix `matrix` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with diagonal Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

## Example
```julia-repl
julia> activity_to_correlation!(C, J=J₀, fourier_type=FastFourier.DCT);
```
"""
function activity_to_correlation!(matrix::AbstractMatrix, J::Function, δ::Real; fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            local val = 2 - exp( -δ*J(q) ) - exp( -δ*J(k) );
            local val /= J(q) + J(k);
            # Take the proper limit for J(k)->0, J(q)->0!
            if !isfinite(val) val = δ end
            matrix[id] *= val;
        end
    end
end

"""
activity_to_correlation!(matrix::AbstractMatrix; J::Function, δ::Real, τ::Real; fourier_type)

This method is used to calculate the pairwise velocity correlation between different points, with lag time τ and measurement window δ. Needs to be divided by δ^2 in the end.

Map the noise correlation matrix `matrix` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with diagonal Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

## Example
```julia-repl
julia> activity_to_correlation!(C, J=J₀, fourier_type=FastFourier.DCT);
```
"""
function activity_to_correlation!(matrix::AbstractMatrix, J::Function, δ::Real, τ::Real; fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            if τ >= δ
                local val = 2exp( -τ*J(q) ) - exp( -(τ+δ)*J(q) ) - exp( -(τ-δ)*J(q) )
                local val /= J(q) + J(k);
                # Take the proper limit for J(k)->0, J(q)->0!
                if !isfinite(val) val = 0 end
            else
                local val = 2exp( -τ*J(q) ) - exp( -(τ+δ)*J(q) ) - exp( -(δ-τ)*J(k) )
                local val /= J(q) + J(k);
                # Take the proper limit for J(k)->0, J(q)->0!
                if !isfinite(val) val = δ-τ end
            end
            matrix[id] *= val;
        end
    end
end

"""
activity_to_correlation!(matrix::AbstractMatrix; J::Function, dJ_dα::Function, fourier_type)

This is a gradient method that computes the response to changing a parameter in the Jacobian. Takes the derivative of the Jacobian as keyword argument.

Map the noise correlation matrix `matrix` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with diagonal Jacobian `J::Function`. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT.

## Example
```julia-repl
julia> activity_to_correlation!(C, J=J₀, fourier_type=FastFourier.DCT);
```
"""
function activity_to_correlation!(matrix::AbstractMatrix, J::Function, dJ_dα::Function; fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            local val = -(dJ_dα(q) + dJ_dα(k));
            local val /= (J(q) + J(k))^2;
            #
            if !isfinite(val) val = 0 end
            matrix[id] *= val
        end
    end
end

"""
correlation_to_activity!(matrix::AbstractMatrix; J::Function, fourier_type)

Map the position correlation in Fourier space `matrix` to return the corresponding noise correlation matrix C of the polymer. Specify `fourier_type` to switch between different fourier transforms such as DCT or FFT. Specify `J` to define the Jacobian of the polymer. 
"""
function correlation_to_activity!(matrix::AbstractMatrix, J::Function; fourier_type)
    for id in CartesianIndices(matrix)
        let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
            # manipulate as needed in Fourier space
            matrix[id] *= J(q) + J(k);
        end
    end
end

end
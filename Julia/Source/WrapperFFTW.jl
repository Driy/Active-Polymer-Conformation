module WrapperFFTW

using FFTW

"""
Specific error for not invalid or not implemented transform types.
"""
struct InvalidTransform <:Exception
    msg::String
end

"""
Wrapper function for forward Fourier transform.
"""
function forward(matrix::AbstractMatrix; fourier_type)
    if fourier_type==DCT
        return dct(matrix)
    elseif fourier_type==FFT
        return fft(matrix)
    else
        InvalidTransform("DCT and FFT are valid fourier_type specifications.")
    end
end

"""
Wrapper function for inplace forward Fourier transform.
"""
function forward!(matrix::AbstractMatrix; fourier_type)
    if fourier_type==DCT
        return dct!(matrix)
    elseif fourier_type==FFT
        return fft!(matrix)
    else
        InvalidTransform("DCT and FFT are valid fourier_type specifications.")
    end
end

"""
Wrapper function for inverse Fourier transform.
"""
function backward(matrix::AbstractMatrix; fourier_type)
    if fourier_type==DCT
        return idct(matrix)
    elseif fourier_type==FFT
        return ifft(matrix)
    else
        InvalidTransform("DCT and FFT are valid fourier_type specifications.")
    end
end

"""
Wrapper function for inplace inverse Fourier transform.
"""
function backward!(matrix::AbstractMatrix; fourier_type)
    if fourier_type==DCT
        return idct!(matrix)
    elseif fourier_type==FFT
        return ifft!(matrix)
    else
        InvalidTransform("DCT and FFT are valid fourier_type specifications.")
    end
end

"""
Define enums that represent the type of Fourier transform.
"""
@enum Type DCT=0 FFT=1

"""
get_frequency(fourier_type, matrix::AbstractMatrix)

Retrieve all Fourier modes for Transform with type `fourier_type`, of the `matrix::AbstractMatrix`.
"""
function frequency(matrix::AbstractMatrix; fourier_type)
    return map(I -> frequency(fourier_type, I, matrix), CartesianIndices(matrix))
end

"""
get_frequency(fourier_type, cartesian_index::CartesianIndex, matrix::AbstractMatrix)

Retrieve Fourier modes for Transform with type `fourier_type`, where `cartesian_index::CartesianIndex` refers to the respective index of the `matrix::AbstractMatrix`.
"""
function frequency(cartesian_index::CartesianIndex, matrix::AbstractMatrix; fourier_type)
    if fourier_type==DCT
        return _frequency_dct(cartesian_index, matrix)
    elseif fourier_type==FFT
        return _frequency_fft(cartesian_index, matrix)
    else
        InvalidTransform("DCT and FFT are valid fourier_type specifications.")
    end
end

"""
get_frequency_dct(cartesian_index::CartesianIndex, matrix::AbstractMatrix)

Retrieve Fourier modes for Discrete Cosine Transform, where `cartesian_index::CartesianIndex` refers to the respective index of the `matrix::AbstractMatrix`.

## Example
```julia-repl
julia> get_frequency_dct(CartesianIndex(1,1), A)
(0.0, 0.0)
```
"""
function _frequency_dct(cartesian_index::CartesianIndex, matrix::AbstractMatrix)
    # this is checked to yield the correct wave modes!
    return π .* (cartesian_index.I .- 1) ./ size(matrix);
end

"""
get_frequency_fft(cartesian_index::CartesianIndex, matrix::AbstractMatrix)

Retrieve Fourier modes for Fast Fourier Transform, where `cartesian_index::CartesianIndex` refers to the respective index of the `matrix::AbstractMatrix`.

## Example
```julia-repl
julia> get_frequency_fft(CartesianIndex(1,1), A)
(0.0, 0.0)
```
"""
function _frequency_fft(cartesian_index::CartesianIndex, matrix::AbstractMatrix)
    # this is checked to yield the correct wave modes!
    return 2π .* (0.5 .- abs.(0.5 .- (cartesian_index.I .- 1) ./ size(matrix) ));
end

end
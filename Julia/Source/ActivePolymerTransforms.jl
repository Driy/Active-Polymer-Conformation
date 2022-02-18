module ActivePolymerTransforms

using FFTW

export @default
export @dct
export @fft
export mean_square_separation
export extract_correlations

"""
J₀(q)

Define Jacobian for a Rouse chain, where `q` is the wave mode. Consider only line tension in Fourier space.

## Example
```julia-repl
julia> J₀(1.0)
1.0
```
"""
function J₀(q)
    return q^2;
end

"""
get_frequency_dct(cartesian_index::CartesianIndex, matrix)

Retrieve Fourier modes for Discrete Cosine Transform, where `cartesian_index::CartesianIndex` refers to the respective index of the `matrix`.

## Example
```julia-repl
julia> get_frequency_dct(CartesianIndex(1,1), A)
(0.0, 0.0)
```
"""
function get_frequency_dct(cartesian_index::CartesianIndex, matrix)
    # this is checked to yield the correct wave modes!
    return π .* (cartesian_index.I .- 1) ./ size(matrix);
end

"""
get_frequency_fft(cartesian_index::CartesianIndex, matrix)

Retrieve Fourier modes for Fast Fourier Transform, where `cartesian_index::CartesianIndex` refers to the respective index of the `matrix`.

## Example
```julia-repl
julia> get_frequency_fft(CartesianIndex(1,1), A)
(0.0, 0.0)
```
"""
function get_frequency_fft(cartesian_index::CartesianIndex, matrix)
    # this is checked to yield the correct wave modes!
    return 2π .* (0.5 .- abs.(0.5 .- (cartesian_index.I .- 1) ./ size(matrix) ));
end

"""
@default

MACRO: Run in default mode. Prepend this to a command.

## Example
```julia-repl
julia> ΔR = @default mean_square_separation(C);
```
"""
macro default(command)
    quote
        $(esc(command))()
    end
end

"""
@dct

MACRO: Run in DCT mode. Prepend this to a command.

## Example
```julia-repl
julia> ΔR = @dct mean_square_separation(C);
```
"""
macro dct(command)
    quote
        $(esc(command))(
            transform_forward = dct, transform_backward = idct, 
            transform_forward! = dct!, transform_backward! = idct!, 
            transform_frequency = get_frequency_dct)
    end
end

"""
@fft

MACRO: Run in FFT mode. Prepend this to a command.

## Example
```julia-repl
julia> ΔR = @fft mean_square_separation(C);
```
"""
macro fft(command)
    quote
        $(esc(command))(
            transform_forward = fft, transform_backward = ifft, 
            transform_forward! = fft!, transform_backward! = ifft!, 
            transform_frequency = get_frequency_fft)
    end
end

"""
fourier_C_to_R!(tmp, n, J::Function, get_frequency::Function)

Map the noise correlation matrix `tmp` in Fourier space to a correlation matrix between different Rouse modes at the same time, for a polymer with Jacobian `J::Function`. Takes the method `get_frequency::Function` to correctly determine the mode numbers of the Fourier modes

## Example
```julia-repl
julia> fourier_C_to_R!(C, 0, J₀, get_frequency_dct);
```
"""
function fourier_C_to_R!(tmp, J::Function, get_frequency::Function)
    for id in CartesianIndices(tmp)
        let (q,k) = get_frequency(id, tmp)
            # manipulate as needed in Fourier space
            tmp[id] /= J(q) + J(k);
        end
    end
    # make sure that the homogeneous mode is not NaN! It is irrelevant for the mean square separation anyways.
    tmp[begin,:] .= 0.0;
    tmp[:,begin] .= 0.0;
end

"""
real_R_to_ΔR(tmp)

Calculate mean square separation matrix from correlation matrix `tmp` between different Rouse modes.

## Example
```julia-repl
julia> ΔR = real_R_to_ΔR(R);
```
"""
function real_R_to_ΔR(tmp)
    return [
        let (x,y) = id.I
            tmp[x,x] + tmp[y,y] - 2tmp[x,y]
        end
        for id in CartesianIndices(tmp)];
end

"""
Calculate the mean square separation of a polymer with given correlation matrix C and Jacobian J.
The keyword arguments specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function mean_square_separation(
        C, J::Function = J₀)
        function( ; # this makes the remaining arguments keyword args!
            transform_forward::Function = dct, transform_backward::Function = idct, 
            transform_forward!::Function = dct!, transform_backward!::Function = idct!, 
            transform_frequency::Function = get_frequency_dct)
        tmp = transform_forward(C);                   # transform to Fourier space; Fourier output is always dense!
        fourier_C_to_R!(tmp, J, transform_frequency); # manipulate as needed in Fourier space
        transform_backward!(tmp);                     # transform to real space
        return real(real_R_to_ΔR(tmp));               # return mean square separation
    end
end

"""
fourier_ΔR_to_C!(tmp, J::Function, get_frequency::Function)

Map the mean square separation in Fourier space `tmp` to `return` a corresponding noise correlation matrix C.
"""
function fourier_ΔR_to_C!(tmp, J::Function, get_frequency::Function)
    for id in CartesianIndices(tmp)
        let (q,k) = get_frequency(id, tmp)
            # manipulate as needed in Fourier space
            tmp[id] *= -0.5;
            tmp[id] *= J(q) + J(k);
        end
    end
    # the homogeneous modes cannot be determined, because ΔR is invariant with respect to them!
    tmp[begin,:] .= 0.0;
    tmp[:,begin] .= 0.0;
end

"""
Extract correlation function from mean square separation `ΔR`, for a polymer with Jacobian `J::Function`.
"""
function extract_correlations(
        ΔR, J::Function = J₀)
    function( ; # this makes the remaining arguments keyword args!
            transform_forward::Function = dct, transform_backward::Function = idct, 
            transform_forward!::Function = dct!, transform_backward!::Function = idct!, 
            transform_frequency::Function = get_frequency_dct)
        tmp = transform_forward(ΔR)                    # transform to Fourier space; Fourier output is always dense!
        fourier_ΔR_to_C!(tmp, J, transform_frequency); # manipulate as needed in Fourier space
        transform_backward!(tmp);                      # transform to real space
        return real(tmp);
    end
end

"""
test_extraction(C)

Test the extraction of the correlation function for a given correlation function.
Note that we cannot extract information about homogeneous contributions!
Thus, this ONLY works perfectly if there are no homogeneous contributions.
"""
function test_extraction(C)
    R = mean_square_separation(C);
    C_extracted = extract_correlations(R);
    return (mean(C - C_extracted), std(C - C_extracted));
end

end
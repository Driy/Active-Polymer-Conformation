module ActivePolymerTransforms

using FFTW

export @default
export @dct
export @fft
export mean_square_separation
export extract_correlations

"""
Define Jacobian for a Rouse chain. Consider only line tension in Fourier space.
"""
function J₀(q)
    return q^2;
end

"""
Retrieve Fourier modes for Discrete Cosine Transform.
"""
function get_frequency_dct(cartesian_index, matrix)
    # this is checked to yield the correct wave modes!
    return π .* (cartesian_index.I .- 1) ./ size(matrix);
end

"""
Retrieve Fourier modes for Fast Fourier Transform.
"""
function get_frequency_fft(cartesian_index, matrix)
    # this is checked to yield the correct wave modes!
    return 2π .* (0.5 .- abs.(0.5 .- (cartesian_index.I .- 1) ./ size(matrix) ));
end

"""
MACRO: Run in default mode.
"""
macro default(command)
    quote
        $(esc(command))()
    end
end

"""
MACRO: Run in DCT mode.
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
MACRO: Run in FFT mode.
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
Map the noise correlation matrix C to a correlation matrix between different Rouse modes at the same time, for a polymer with Jacobian J.
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
Calculate mean square separation matrix from correlation matrix between different Rouse modes.
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
Map the mean square separation in Fourier space to a corresponding noise correlation matrix C.
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
Extract correlation function from mean square separation, for a polymer with Jacobian J.
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
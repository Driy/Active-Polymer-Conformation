module CorrelationMatrices

using FFTW

export random

"""
Generate random correlation matrix of size `matrix_size`. 
If `remove_homogeneous` is true, then remove homogeneous contributions.
"""
function random(matrix_size; remove_homogeneous=false)
    # first we define the corresponding random matrix Λ
    tmp = rand(matrix_size, matrix_size) .- 0.5;
    
    # remove homogeneous contributions if requested
    if remove_homogeneous
        dct!(tmp); tmp[begin,:] .= 0.0; tmp[:, begin] .= 0.0; idct!(tmp);
    end

    # return positive definite correlation matrix C
    return tmp*tmp';
end

end
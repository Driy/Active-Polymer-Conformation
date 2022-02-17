module CorrelationMatrices

using FFTW

export random_correlation_matrix
export random_correlation_matrix_special

"""
Generate random correlation matrix of given size.
"""
function random_correlation_matrix(matrix_size)
    # first we define the corresponding random matrix Λ
    tmp = rand(matrix_size, matrix_size) .- 0.5;

    # return positive definite correlation matrix C
    return tmp*tmp';
end

"""
Generate random correlation matrix of given size, which has no homogeneous contributions.
"""
function random_correlation_matrix_special(matrix_size)
    # first we define the corresponding random matrix Λ
    tmp = rand(matrix_size, matrix_size) .- 0.5;

    # remove homogeneous contributions
    dct!(tmp); tmp[begin,:] .= 0.0; tmp[:, begin] .= 0.0; idct!(tmp);

    # return positive definite correlation matrix C
    return tmp*tmp';
end

end
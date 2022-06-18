module CorrelationMatrices

using FFTW
using LinearAlgebra

export diagonal_generic, dense_random, dense_box

"""
Generate diagonal correlation matrix of size `matrix_size`. Pass `profile` to set the spatial profile of the independent activity modulations.
"""
function diagonal_generic(matrix_size::Int64; profile::Function=(s->1.0))
    return [profile(s) for s=0:matrix_size-1] |> diagm;
end

"""
Generate random correlation matrix of size `matrix_size`. 
If `remove_homogeneous` is true, then remove homogeneous contributions.
"""
function dense_random(matrix_size::Int64; remove_homogeneous=false)
    # first we define the corresponding random matrix Λ
    tmp = rand(matrix_size, matrix_size) .- 0.5;
    
    # remove homogeneous contributions if requested
    if remove_homogeneous
        dct!(tmp); tmp[begin,:] .= 0.0; tmp[:, begin] .= 0.0; idct!(tmp);
    end

    # return positive definite correlation matrix C
    return tmp*tmp';
end

"""
Generate box correlation matrix of size `matrix_size`.
"""
function dense_box(matrix_size::Int64; box_fraction=0.1)
    vec = [abs(i/matrix_size-0.5)<box_fraction for i in 0:matrix_size-1];    
    return vec*vec';
end

end
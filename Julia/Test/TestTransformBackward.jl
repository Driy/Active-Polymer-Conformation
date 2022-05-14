module TestTransformBackward

include("../Source/ActivePolymer.jl")
include("../Source/CorrelationMatrices.jl")

using Statistics
using LinearAlgebra

"""
test()

Check if inverse transform yields the input of the forward transform for a constant activity. Keyword arguments: `N` specifies the size of the matrix, `tol` the tolerance of the test.
"""
function test_homogeneous(; N=512, tol=1e-2)
    C_diagonal    = CorrelationMatrices.diagonal_generic(N, profile=(s->1.0));
    R_numeric, ΔR_numeric = ActivePolymer.TransformForward.compute_conformation(C_diagonal);
    C_extract     = ActivePolymer.TransformBackward.extract_excitations(R_numeric, ΔR_numeric);

    return abs.(C_diagonal - C_extract) |> mean |> x->(x<tol)
end

"""
test(ϵ, α)

Check if inverse transform yields the input of the forward transform for a constant activity that is superimposed by sinusoidal modulations `ϵ cos(α s)`. Keyword arguments: `N` specifies the size of the matrix, `tol` the tolerance of the test.
"""
function test_cosine(ϵ, α; N=512, tol=1e-2)
    C_diagonal    = CorrelationMatrices.diagonal_generic(N, profile=(s->1.0 + ϵ*cos(α*s)));
    R_numeric, ΔR_numeric = ActivePolymer.TransformForward.compute_conformation(C_diagonal);
    C_extract     = ActivePolymer.TransformBackward.extract_excitations(R_numeric, ΔR_numeric);

    return abs.(C_diagonal - C_extract) |> mean |> x->(x<tol)
end

"""
test(ϵ, α)

Check if inverse transform yields the input of the forward transform for a random activity matrix. Keyword arguments: `N` specifies the size of the matrix, `tol` the tolerance of the test, `num` the number of tests, `remove_homogeneous` if the random activity matrix is devoid of homogeneous modes.
"""
function test_random(; N=512, tol=1e-2, num=1, remove_homogeneous=false)
    function inner()
        C_random      = CorrelationMatrices.dense_random(N, remove_homogeneous = remove_homogeneous);
        R_numeric, ΔR_numeric = ActivePolymer.TransformForward.compute_conformation(C_random);
        C_extract     = ActivePolymer.TransformBackward.extract_excitations(R_numeric, ΔR_numeric);
        return abs.((C_random - C_extract) / Diagonal(C_random) ) |> mean |> x->(x<tol)
    end

    return all([inner() for i=1:num])
end

end
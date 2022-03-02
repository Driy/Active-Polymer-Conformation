module TestForwardTransform

include("../Source/ActivePolymer.jl")
include("../Source/CorrelationMatrices.jl")

using Statistics
using LinearAlgebra

"""
test()

Check if discrete computation yields the same result as the analytic prediction for a constant activity. Keyword arguments: `N` specifies the size of the matrix, `tol` the tolerance of the test.
"""
function test_homogeneous(; N=512, tol=1e-2)
    C_diagonal  = CorrelationMatrices.diagonal_generic(N, profile=(s->1.0));
    ΔR_numeric  = ActivePolymer.ForwardTransform.mean_square_separation(C_diagonal);
    ΔR_analytic = [analytic(s₁, s₂) for s₁=0:N-1, s₂=0:N-1];

    return abs.(ΔR_numeric - ΔR_analytic) |> mean |> x->(x/N) |> x->(x<tol)
end

"""
test(ϵ, α)

Check if discrete computation yields the same result as the analytic prediction for a constant activity that is superimposed by sinusoidal modulations `ϵ cos(α s)`. Keyword arguments: `N` specifies the size of the matrix, `tol` the tolerance of the test.
"""
function test_cosine_1(ϵ, α; N=512, tol=1e-2)
    C_diagonal  = CorrelationMatrices.diagonal_generic(N, profile=(s->1.0 + ϵ*cos(α*s)));
    ΔR_numeric  = ActivePolymer.ForwardTransform.mean_square_separation(C_diagonal);
    ΔR_analytic = [analytic(ϵ, α, s₁, s₂) for s₁=0:N-1, s₂=0:N-1];

    return abs.(ΔR_numeric - ΔR_analytic) |> mean |> x->(x/N) |> x->(x<tol)
end

"""
Check if discrete computation yields the same result as the analytic prediction for a constant activity that is superimposed by sinusoidal modulations `ϵ cos(α s)`. Keyword arguments: `N` specifies the size of the matrix, `dN` the width of the boundary padding that is to be discarded, `tol` the tolerance of the test.
"""
function test_cosine_2(ϵ, α; N=512, dN=64, tol=1e-2)
    C_diagonal  = CorrelationMatrices.diagonal_generic(N, profile=(s->ϵ*cos(α*s)));
    ΔR_numeric  = ActivePolymer.ForwardTransform.mean_square_separation(C_diagonal);
    ΔR_analytic = [analytic(ϵ, α, s₁, s₂) - analytic(s₁, s₂) for s₁=0:N-1, s₂=0:N-1];

    return abs.(ΔR_numeric - ΔR_analytic)[dN:end-dN,dN:end-dN] |> mean |> x->(α*x) |> x->(x<tol)
end

"""
analytic(s₁, s₂)

Determine the analyic mean squared separation between position s₁ and s₂. The analytic result for an open polymer with constant activity `1` is:
`Δs/2`, where `Δs:=|s₁-s₂|`.
"""
function analytic(s₁, s₂)
    Δs = abs(s₁-s₂);
    return 0.5Δs;
end

"""
analytic(ϵ, α, s₁, s₂)

Determine the analyic mean squared separation between position s₁ and s₂. The analytic result for an open polymer with activity modulations `1 + ϵ cos(α s)` is:
`Δs/2 + (ϵ/α) cos(α S/2) [cos(α Δs/2) - exp(-α Δs/2)]`, where `Δs:=|s₁-s₂|` and `S:=s₁+s₂`.
"""
function analytic(ϵ, α, s₁, s₂)
    Δs = abs(s₁-s₂);
    S  = s₁+s₂;
    return 0.5Δs + (ϵ/α) * cos(α * 0.5S) * ( cos(α * 0.5Δs) - exp(-α * 0.5Δs) )
end

end
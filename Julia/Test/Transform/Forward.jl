module Forward

include("../../Source/ActivePolymer.jl")
include("../../Source/CorrelationMatrices.jl")

using Statistics
using LinearAlgebra

"""
test()

Check if discrete computation yields the same result as the analytic prediction for a constant activity. Keyword arguments: `N` specifies the size of the matrix, `tol` the tolerance of the test.
"""
function test_homogeneous(; N=512, tol=1e-2)
    C_diagonal  = CorrelationMatrices.diagonal_generic(N, profile=(s->1.0));
    ΔR_numeric  = ActivePolymer.Transform.Forward.compute_conformation(C_diagonal) |> ActivePolymer.Methods.Real.correlation_to_separation;
    ΔR_analytic = [ActivePolymer.Methods.Analytic.separation_rouse(s₁, s₂) for s₁=0:N-1, s₂=0:N-1];

    return abs.(ΔR_numeric - ΔR_analytic) |> mean |> x->(x/N) |> x->(x<tol)
end

"""
test(ϵ, α)

Check if discrete computation yields the same result as the analytic prediction for a constant activity that is superimposed by sinusoidal modulations `ϵ cos(α s)`. Keyword arguments: `N` specifies the size of the matrix, `tol` the tolerance of the test.
"""
function test_cosine_1(ϵ, α; N=512, tol=1e-2)
    C_diagonal  = CorrelationMatrices.diagonal_generic(N, profile=(s->1.0 + ϵ*cos(α*s)));
    ΔR_numeric  = ActivePolymer.Transform.Forward.compute_conformation(C_diagonal) |> ActivePolymer.Methods.Real.correlation_to_separation;
    ΔR_analytic = [ActivePolymer.Methods.Analytic.separation_rouse(s₁, s₂, amplitude=ϵ, wavemode=α) for s₁=0:N-1, s₂=0:N-1];

    return abs.(ΔR_numeric - ΔR_analytic) |> mean |> x->(x/N) |> x->(x<tol)
end

"""
Check if discrete computation yields the same result as the analytic prediction for a constant activity that is superimposed by sinusoidal modulations `ϵ cos(α s)`. Keyword arguments: `N` specifies the size of the matrix, `dN` the width of the boundary padding that is to be discarded, `tol` the tolerance of the test.
"""
function test_cosine_2(ϵ, α; N=512, dN=64, tol=1e-2)
    C_diagonal  = CorrelationMatrices.diagonal_generic(N, profile=(s->ϵ*cos(α*s)));
    ΔR_numeric  = ActivePolymer.Transform.Forward.compute_conformation(C_diagonal) |> ActivePolymer.Methods.Real.correlation_to_separation;
    ΔR_analytic = [ActivePolymer.Methods.Analytic.separation_rouse(s₁, s₂, amplitude=ϵ, wavemode=α) - ActivePolymer.Methods.Analytic.separation_rouse(s₁, s₂) for s₁=0:N-1, s₂=0:N-1];

    return abs.(ΔR_numeric - ΔR_analytic)[begin+dN:end-dN,begin+dN:end-dN] |> mean |> x->(α*x) |> x->(x<tol)
end

end
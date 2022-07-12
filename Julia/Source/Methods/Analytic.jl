module Analytic

export separation_rouse, separation_generic

"""
separation_rouse(s₁, s₂; amplitude=0, wavemode=0)

Determine the analyic mean squared separation between the Rouse monomers at position s₁ and s₂. The analytic result for an open Rouse polymer with activity modulations `1 + ϵ cos(α s)` is:
`Δs/2 + (ϵ/α) cos(α S/2) [cos(α Δs/2) - exp(-α Δs/2)]`, where `Δs:=|s₁-s₂|` and `S:=s₁+s₂`.
"""
function separation_rouse(s₁, s₂; amplitude::Real=0, wavemode::Real=1)
    Δs = abs(s₁-s₂);
    S  = s₁+s₂;
    return 0.5Δs + (amplitude/wavemode) * cos(wavemode * 0.5S) * ( cos(wavemode * 0.5Δs) - exp(-wavemode * 0.5Δs) )
end

"""
correlation_generic(s₁, s₂; amplitude=0, wavemode=0)

Determine the analyic correlation between the monomers at position s₁ and s₂.
"""
function correlation_generic(Δs; T::Real, α::Real, κ::Real)
    # compute analytic solution for homogeneous activity matrix
    determinant = 1 - 4α * κ |> ComplexF64 |> sqrt
    tmp = sqrt(1+determinant)*exp.(-Δs*sqrt( (1-determinant)/(2κ) ))
    tmp -= sqrt(1-determinant)*exp.(-Δs*sqrt( (1+determinant)/(2κ) ))
    tmp /= 4sqrt(2α)*determinant
    return T*tmp |> real
end

"""
separation_generic(s₁, s₂; amplitude=0, wavemode=0)

Determine the analyic mean squared separation between the monomers at position s₁ and s₂.
"""
function separation_generic(Δs; parameters...)
    return 2( correlation_generic(0.0; parameters...) .- correlation_generic(Δs; parameters...) )
end

end
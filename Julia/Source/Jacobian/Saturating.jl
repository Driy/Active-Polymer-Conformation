module Saturating

"""
J(q; α=0, κ=0, n=4)

Define Jacobian for a Rouse chain, where `q` is the wave mode. Consider only line tension in Fourier space.
Measured in terms of characteristic activity over characteristic length squared times inverse friction.

## Example
```julia-repl
julia> J(1.0)
1.0
```
"""
function J(α::Real=0.0, κ::Real=0.0, n::Real=2)
    @. return function(q)
        return 0.5 * (α + q^2) / (1.0 + κ * abs(q)^n);
    end
end

"""
dJ_dα(q; α=0, κ=0, n=4)

Define Jacobian for a Rouse chain, where `q` is the wave mode. Consider only line tension in Fourier space.
Measured in terms of characteristic activity over characteristic length squared times inverse friction.

## Example
```julia-repl
julia> dJ_dα(1.0)
1.0
```
"""
function dJ_dα(α::Real=0.0, κ::Real=0.0, n::Real=2)
    @. return function(q)
        return (0.5) / (1.0 + κ * abs(q)^n);
    end
end

"""
dJ_dκ(q; α=0, κ=0, n=4)

Define Jacobian for a Rouse chain, where `q` is the wave mode. Consider only line tension in Fourier space.
Measured in terms of characteristic activity over characteristic length squared times inverse friction.

## Example
```julia-repl
julia> dJ_dκ(1.0)
1.0
```
"""
function dJ_dκ(α::Real=0.0, κ::Real=0.0, n::Real=2)
    @. return function(q)
        return -0.5*(α + q^2) / (1.0 + κ * abs(q)^n)^2 * abs(q)^n;
    end
end

end
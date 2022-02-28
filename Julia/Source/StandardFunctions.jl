module StandardFunctions

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

end
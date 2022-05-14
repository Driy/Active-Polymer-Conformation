module MethodsReal

using LinearAlgebra

using ..WrapperFFTW

export correlation_to_position, correlation_to_separation, position_and_separation_to_correlation

"""
correlation_to_position(matrix::AbstractMatrix)

Calculate the average radial position from the correlation matrix between different Rouse modes (`matrix`).

## Example
```julia-repl
julia> ΔR = X_to_ΔR(R);
```
"""
function correlation_to_position(matrix::AbstractMatrix)
    return matrix |> diag |> real;
end

"""
correlation_to_separation(matrix::AbstractMatrix)

Calculate mean squared separation matrix from the correlation matrix between different Rouse modes (`matrix`).

## Example
```julia-repl
julia> ΔR = X_to_ΔR(R);
```
"""
function correlation_to_separation(matrix::AbstractMatrix)
    return [
        let (x,y) = id.I
            matrix[x,x] + matrix[y,y] - 2matrix[x,y]
        end
        for id in CartesianIndices(matrix)] |> real;
end

"""
position_and_separation_to_correlation(matrix::AbstractMatrix)

Calculate mean squared separation matrix from the correlation matrix between different Rouse modes (`matrix`).

## Example
```julia-repl
julia> ΔR = real_R_to_ΔR(R);
```
"""
function position_and_separation_to_correlation(vector_position::AbstractVector, matrix_separation::AbstractMatrix)
    return [
        let (x,y) = id.I
            0.5( vector_position[x] + vector_position[y] - matrix_separation[x,y] )
        end
        for id in CartesianIndices(matrix_separation)] |> real;
end

end
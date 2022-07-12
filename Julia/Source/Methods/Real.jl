module Real

using LinearAlgebra
using Statistics

export correlation_to_position, correlation_to_separation, position_and_separation_to_correlation, correlation_split, marginalize_translation

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
julia> ΔR = correlation_to_separation(RX);
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
correlation_split(matrix::AbstractMatrix)

Calculate average radial position and mean squared separation matrix from the correlation matrix between different Rouse modes (`matrix`).

## Example
```julia-repl
julia> ΔR = correlation_to_separation(RX);
```
"""
function correlation_split(matrix::AbstractMatrix)
    return correlation_to_position(matrix), correlation_to_separation(matrix)
end

"""
position_and_separation_to_correlation(matrix::AbstractMatrix)

Calculate mean squared separation matrix from the correlation matrix between different Rouse modes (`matrix`).

## Example
```julia-repl
julia> RX = position_and_separation_to_correlation(R, ΔR);
```
"""
function position_and_separation_to_correlation(vector_position::AbstractVector, matrix_separation::AbstractMatrix)
    return [
        let (x,y) = id.I
            0.5( vector_position[x] + vector_position[y] - matrix_separation[x,y] )
        end
        for id in CartesianIndices(matrix_separation)] |> real;
end

"""
marginalize_translation(matrix::AbstractMatrix)

Marginalize diagonal translations of the given matrix. Returns average value as a function of the distance to the diagonal.
"""
function marginalize_translation(matrix::AbstractMatrix)
    return [diag(matrix,i) |> mean for i in 0:size(matrix,1)-1];
end

end
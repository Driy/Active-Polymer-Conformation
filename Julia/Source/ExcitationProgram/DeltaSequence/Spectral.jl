module Spectral

"""
Takes the delay sequence and the excitation sequence in Fourier space, and gives the effective resulting correlation (that is, effective correlations between instantaneous excitations that would yield the same correlation between Rouse modes).
"""
function effective_correlation(delay_sequence, excitation_sequence, J::Function; fourier_type)
    @assert length(delay_sequence) == length(excitation_sequence)
    
    # initialize matrix which will be returned
    matrix = zeros(excitation_sequence[begin] |> size)
    
    # iterate over all excitations
    for α, C in delay_sequence, excitation_sequence
        # iterate over all matrix elements
        for id in CartesianIndices(matrix)
            let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
                # manipulate as needed in Fourier space
                tmp = exp( -α*J(q) ) * C'[id] + C[id] * exp( -α*J(k) )
                matrix[id] += tmp                
            end
        end        
    end
    return matrix
end

"""
Takes the delay sequence and the excitation sequence in Fourier space, and gives the correction term ('B') for the fluctuations.
"""
function delayed_correlation(Δt, delay_sequence, excitation_sequence, J::Function; fourier_type)
    @assert length(delay_sequence) == length(excitation_sequence)
    
    # initialize matrix which will be returned
    matrix = zeros(excitation_sequence[begin] |> size)
    
    # iterate over all excitations
    for α, C in delay_sequence, excitation_sequence
        # iterate over all matrix elements
        for id in CartesianIndices(matrix)
            let (q,k) = FastFourier.frequency(id, matrix, fourier_type = fourier_type)
                # manipulate as needed in Fourier space
                tmp = C[id] * exp( -α*J(k) )
                tmp /= J(q) + J(k)
                tmp *= exp( (J(q) + J(k)) * min(Δt, α) ) - 1.0
                
                # check if we divided by zero and, if yes, then take the proper limit
                if !isfinite(tmp)
                    tmp = C[id] * min(Δt, α)
                end
                matrix[id] += tmp                
            end
        end        
    end
    return matrix
end

end
module Real

using LinearAlgebra

using ..Spectral
using ....Jacobian
using ....Methods
using ....Methods.FastFourier

"""
compute_pairwise_velocity_correlation(matrix_activity::AbstractMatrix; J::Function = Jacobian.Standard.J₀, fourier_type = FastFourier.DCT)

Calculate the mean squared traveled distance after the time τ, for all loci along a polymer with given correlation matrix `matrix_activity` and Jacobian `J`.
Use the argument `fourier_type` to specialize this function to specific transforms, such as the Discrete Cosine Transform or the Fast Fourier Transform.
"""
function compute_pairwise_velocity_correlation(delay_sequence, excitation_sequence, J::Function, δ::Real, τ::Real;
        fourier_type = FastFourier.DCT)
    # transform to Fourier space; always dense!
    excitation_sequence_fourier = [FastFourier.forward(matrix_activity, fourier_type = fourier_type) 
        for matrix_activity in excitation_sequence];
    
    # get effective excitation matrix
    tmp = Spectral.effective_correlation(delay_sequence, excitation_sequence_fourier, J; fourier_type = fourier_type)

    # manipulate as needed in Fourier space
    Methods.Spectral.activity_to_correlation!(tmp, J, δ, τ, fourier_type = fourier_type);
    Spectral.apply_correlation_correction!(
        tmp, delay_sequence, excitation_sequence_fourier, J, δ, τ, fourier_type = fourier_type);
    
    # transform to real space
    FastFourier.backward!(tmp, fourier_type = fourier_type);
    # return mean square separation
    return tmp ./ δ^2;
end

end
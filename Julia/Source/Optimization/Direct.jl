module Direct

using HDF5
using LinearAlgebra
using Statistics
using Distributed
using SharedArrays
using ProgressMeter

using ...Methods
using ...Jacobian
using ...CorrelationMatrices

using ..Model
using ..Interface

"""
coupling_matrix!(mat, i, j, parameters; kwargs...)

Determine elements of the coupling matrix, which determines the gradient of the [proposed mean squared separation matrix] squared, with respect to localized changes in activity.
"""
function populate_coupling_matrix!(mat, i, j, parameters; kwargs...)
    N = size(mat, 1)
    tmp1 = Model.numeric_dense(
        CorrelationMatrices.diagonal_delta(i, N), parameters[2:end]...; kwargs...);
    tmp2 = Model.numeric_dense(
        CorrelationMatrices.diagonal_delta(j, N), parameters[2:end]...; kwargs...);
    mat[i,j] = mean(tmp1 .* tmp2);
end

"""
coupling_vector!(mat, i, j, parameters; kwargs...)

Determine elements of the coupling vector, which determines the gradient of the [proposed mean squared separation matrix] times the [mean squared separation data], with respect to localized changes in activity.
"""
function populate_coupling_vector!(vec, i, ΔR, parameters; kwargs...)
    N = size(vec, 1)
    tmp1 = Model.numeric_dense(
        CorrelationMatrices.diagonal_delta(i, N), parameters[2:end]...; kwargs...);
    vec[i] = mean(tmp1 .* ΔR);
end

function setup_direct_system(name; jacmodule=Jacobian.Discrete, modeltype=Model.Full, n=3)
    R, ΔR = Interface.load_data(name);
    RX = Methods.Real.position_and_separation_to_correlation(R, ΔR);
    N = size(ΔR, 1)

    # Three-parameter fits
    parameters  = Interface.fit_mechanics(
        ΔR, ModelType=modeltype, jacmodule=jacmodule, n=n)

    # get rhs vector v in linear equation
    vec = SharedVector{Float64}(N)
    @showprogress @distributed for i in 1:N
        populate_coupling_vector!(vec, i, ΔR, parameters.minimizer, jacmodule=jacmodule, n=n)
    end

    # get lhs matrix M in linear equation. Solution (optimal distribution of activity) is M^-1 * v
    mat = SharedMatrix{Float64}(N,N)
    @showprogress @distributed for i in 1:N
        for j in 1:i
            populate_coupling_matrix!(mat, i, j, parameters.minimizer, jacmodule=jacmodule, n=n)
        end
    end
    matsym = Symmetric(mat, :L);
    
    if jacmodule==Jacobian.Discrete
        jacobian_type="discrete"
    elseif jacmodule==Jacobian.Standard
        jacobian_type="standard"
    end

    file = h5open(["data/", name, "_", jacobian_type, "_n=", n , ".hdf5"] |> join, "w")
    create_group(file, "optimization")
    file["optimization/matrix"]=Matrix(matsym)
    file["optimization/vector"]=vec
    close(file)
end

end
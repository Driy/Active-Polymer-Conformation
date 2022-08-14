module DirectComparable

using HDF5
using Statistics
using Distributed
using SharedArrays
using LinearAlgebra
using ProgressMeter

using ...Methods
using ...Jacobian
using ...CorrelationMatrices

using ..Model
using ..Interface

"""
coupling_matrix(N, parameters; kwargs...)

Determine elements of the coupling matrix, which determines the gradient of the [proposed mean squared separation matrix] squared, with respect to localized changes in activity.
"""
function coupling_matrix(N, parameters; kwargs...)
    mat = SharedMatrix{Float64}(N,N);
    @showprogress @distributed for i in 0:N-1
        for j in 0:min(i, N-1-i)
            # calculate matrix elements
            tmp1 = Model.numeric_dense(
                CorrelationMatrices.diagonal_delta(1+i, N), parameters[2:end]...; kwargs...);
            tmp2 = Model.numeric_dense(
                CorrelationMatrices.diagonal_delta(1+j, N), parameters[2:end]...; kwargs...);
            val = mean(tmp1 .* tmp2);

            # populate matrix and exploit its fourfold symmetry
            mat[end-i,end-j]=val;
            mat[end-j,end-i]=val;
            mat[begin+i,begin+j]=val;
            mat[begin+j,begin+i]=val;
        end
        # sleep necessary for progress bar update
        sleep(1e-4)
    end
    return mat;
end

"""
coupling_vector!(mat, i, j, parameters; kwargs...)

Determine elements of the coupling vector, which determines the gradient of the [proposed mean squared separation matrix] times the [mean squared separation data], with respect to localized changes in activity.
"""
function coupling_vector(ΔR, parameters; kwargs...)
    N = size(ΔR, 1)
    vec = SharedVector{Float64}(N);
    @showprogress @distributed for i in 0:N-1
        # calculate vector elements
        tmp1 = Model.numeric_dense(
            CorrelationMatrices.diagonal_delta(1+i, N), parameters[2:end]...; kwargs...);
        val = mean(tmp1 .* ΔR);
        
        # populate vector
        vec[begin+i] = val;
        
        # sleep necessary for progress bar
        sleep(1e-4)
    end
    return vec;
end

function setup_reference_system(name; jacmodule=Jacobian.Discrete, modeltype=Model.Full, n=3, padding::Real=0.85)
    
    # load data
    R, ΔR = Interface.load_data(name);

    # Three-parameter fits
    parameters  = Interface.fit_mechanics(
        ΔR, modeltype=modeltype, jacmodule=jacmodule, n=n, padding=padding)
    
    # get lhs matrix M in linear equation. Solution (optimal distribution of activity) is M^-1 * v
    mat = coupling_matrix(size(ΔR, 1), parameters.minimizer, jacmodule=jacmodule, n=n);
    
    return parameters, mat, jacmodule, modeltype, n
end

function setup_direct_system(name, reference; overwrite=false)
    
    # import reference system
    parameters, mat, jacmodule, modeltype, n = reference
    
    # 
    if jacmodule==Jacobian.Discrete
        jacobian_type="discrete"
    elseif jacmodule==Jacobian.Standard
        jacobian_type="standard"
    end
    
    # check if file exists and decide what to do
    if (["data/", name, "_", jacobian_type, "_n=", n , ".hdf5"] |> join |> isfile) && overwrite==false
        return
    end
    
    R, ΔR = Interface.load_data(name);

    # get rhs vector v in linear equation
    vec = coupling_vector(ΔR, parameters.minimizer, jacmodule=jacmodule, n=n);
        
    file = h5open(["data/", name, "_", jacobian_type, "_n=", n , ".hdf5"] |> join, "w")
    create_group(file, "optimization")
    file["optimization/matrix"]=Matrix(mat)
    file["optimization/vector"]=Vector(vec)
    close(file)
end

end
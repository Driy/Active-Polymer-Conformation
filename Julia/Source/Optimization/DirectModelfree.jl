module DirectModelfree

using HDF5
using FFTW
using Statistics
using Distributed
using SharedArrays
using LinearAlgebra
using ProgressMeter
using ImageFiltering

using ...Methods
using ...Transform
using ...CorrelationMatrices

using ..Model
using ..Interface

"""
coupling_matrix(N, parameters; kwargs...)

Determine elements of the coupling matrix, which determines the gradient of the [proposed mean squared separation matrix] squared, with respect to localized changes in activity.
"""
function coupling_matrix(J)
    N = size(J, 1)
    mat = SharedMatrix{Float64}(N,N);
    @showprogress @distributed for i in 0:N-1
        for j in 0:min(i, N-1-i)
            # calculate matrix elements
            
            tmp1 = Transform.Forward.compute_conformation(
                CorrelationMatrices.diagonal_delta(1+i, N), J) |> Methods.Real.correlation_to_separation;
            tmp2 = Transform.Forward.compute_conformation(
                CorrelationMatrices.diagonal_delta(1+j, N), J) |> Methods.Real.correlation_to_separation;            
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
function coupling_vector(ΔR, J)
    N = size(ΔR, 1)
    vec = SharedVector{Float64}(N);
    @showprogress @distributed for i in 0:N-1
        # calculate vector elements
        tmp1 = Transform.Forward.compute_conformation(
            CorrelationMatrices.diagonal_delta(1+i, N), J) |> Methods.Real.correlation_to_separation;        
        val = mean(tmp1 .* ΔR);
        
        # populate vector
        vec[begin+i] = val;
        
        # sleep necessary for progress bar
        sleep(1e-4)
    end
    return vec;
end

function setup_direct_system(name; overwrite=false, window=2)
        
    # check if file exists and decide what to do
    if (["data/", name, "_modelfree", ".hdf5"] |> join |> isfile) && overwrite==false
        return
    end
    
    R, ΔR = Interface.load_data(name);

    # extract model
    ΔR_marginalized = Methods.Real.marginalize_translation(ΔR);
    N = size(ΔR_marginalized,1)
    J = -1 ./ (sqrt(N) * dct(ΔR_marginalized)) |> x->imfilter(x, Kernel.gaussian((window,)))    

    # get rhs vector v in linear equation
    vec = coupling_vector(ΔR, J);
    
    # get lhs matrix M in linear equation. Solution (optimal distribution of activity) is M^-1 * v
    mat = coupling_matrix(J);
    
    file = h5open(["data/", name, "_modelfree", ".hdf5"] |> join, "w")
    create_group(file, "optimization")
    file["optimization/matrix"]=Matrix(mat)
    file["optimization/vector"]=Vector(vec)
    close(file)
end

end
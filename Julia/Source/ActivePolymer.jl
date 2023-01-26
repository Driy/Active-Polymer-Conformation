module ActivePolymer

include("./CorrelationMatrices.jl")

module Jacobian
include("./Jacobian/Standard.jl")
include("./Jacobian/Discrete.jl")
include("./Jacobian/Saturating.jl")
end

module Methods
include("./Methods/FastFourier.jl")
include("./Methods/Real.jl")
include("./Methods/Spectral.jl")
include("./Methods/Analytic.jl")
end

module Transform
include("./Transform/Forward.jl")
include("./Transform/Backward.jl")
end

module Optimization
include("./Optimization/Model.jl")
include("./Optimization/Residual.jl")
include("./Optimization/Interface.jl")
include("./Optimization/Direct.jl")
include("./Optimization/DirectModelfree.jl")
include("./Optimization/DirectComparable.jl")
end

module ExcitationProgram
module DeltaSequence
include("./ExcitationProgram/DeltaSequence/Analytic.jl")
include("./ExcitationProgram/DeltaSequence/Spectral.jl")
include("./ExcitationProgram/DeltaSequence/RealSpace.jl")
end
end

end
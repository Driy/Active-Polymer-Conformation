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
end

end
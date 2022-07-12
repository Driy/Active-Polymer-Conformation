module ActivePolymer

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

include("./CorrelationMatrices.jl")
include("./MethodsFitting.jl")

end
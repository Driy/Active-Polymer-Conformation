module ActivePolymer

include("./CorrelationMatrices.jl")
include("./JacobianStandard.jl")
include("./JacobianDiscrete.jl")
include("./JacobianSaturating.jl")
include("./WrapperFFTW.jl")
include("./MethodsReal.jl")
include("./MethodsSpectral.jl")
include("./MethodsAnalytic.jl")
include("./TransformForward.jl")
include("./TransformBackward.jl")
include("./MethodsFitting.jl")

export WrapperFFTW, MethodsReal, MethodsSpectral, MethodsAnalytic, TransformForward, TransformBackward, MethodsFitting

end
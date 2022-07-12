module TestTransform
include("./Transform/Forward.jl")
include("./Transform/Backward.jl")
end

using Test

@testset "Validation of transforms." verbose=true begin
    @testset "Forward transforms vs analytic results." begin
        @test TestTransform.Forward.test_homogeneous()
        @test TestTransform.Forward.test_cosine_1(0.5, 12.0pi/512.0)
        @test TestTransform.Forward.test_cosine_2(0.5, 12.0pi/512.0)
    end

    @testset "Inverse transforms vs input values." begin
        @test TestTransform.Backward.test_homogeneous()
        @test TestTransform.Backward.test_cosine(0.5, 12.0pi/512.0)
        @test TestTransform.Backward.test_random(num=64, remove_homogeneous=false)
        @test TestTransform.Backward.test_random(num=64, remove_homogeneous=true)
    end
end

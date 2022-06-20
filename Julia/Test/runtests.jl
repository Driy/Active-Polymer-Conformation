using Test

include("./TestTransformForward.jl")
include("./TestTransformBackward.jl")

@testset "Validation of transforms." verbose=true begin
    @testset "Forward transforms vs analytic results." begin
        @test TestTransformForward.test_homogeneous()
        @test TestTransformForward.test_cosine_1(0.5, 12.0pi/512.0)
        @test TestTransformForward.test_cosine_2(0.5, 12.0pi/512.0)
    end

    @testset "Inverse transforms vs input values." begin
        @test TestTransformBackward.test_homogeneous()
        @test TestTransformBackward.test_cosine(0.5, 12.0pi/512.0)
        @test TestTransformBackward.test_random(num=64, remove_homogeneous=false)
        @test TestTransformBackward.test_random(num=64, remove_homogeneous=true)
    end
end

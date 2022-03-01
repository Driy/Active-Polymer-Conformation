using Test

include("./Test/TestForwardTransform.jl")
include("./Test/TestInverseTransform.jl")

@testset "Validation of transforms." verbose=true begin
    @testset "Forward transforms vs analytic results." begin
        @test TestForwardTransform.test_homogeneous()
        @test TestForwardTransform.test_cosine_1(0.5, 12.0pi/512.0)
        @test TestForwardTransform.test_cosine_2(0.5, 12.0pi/512.0)
    end

    @testset "Inverse transforms vs input values." begin
        @test TestInverseTransform.test_homogeneous()
        @test TestInverseTransform.test_cosine(0.5, 12.0pi/512.0)
        @test TestInverseTransform.test_random(num=64, remove_homogeneous=false)
        @test TestInverseTransform.test_random(num=64, remove_homogeneous=true)
    end
end
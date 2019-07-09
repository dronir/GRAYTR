@testset "Spectra" begin


@testset "NoLight" begin

    N = NoLight()
    
    @test N == GRAYTR.nolight
    @test 2.0 .* N == N
    @test N .* 2.0 == N
    @test N .* N == N
    @test N .+ N == N
    

end



@testset "SampledSpectrum" begin

    S1 = SampledSpectrum(300, 800, [0.0, 1.0, 0.0])
    S2 = SampledSpectrum(300, 800, [0.5, 0.5, 0.5])
    
    @test S1.values ≈ [0.0, 1.0, 0.0]
    @test S2.values ≈ [0.5, 0.5, 0.5]

    # SampledSpectrum with itself
    
    S3 = S1 .+ S2
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [0.5, 1.5, 0.5]
    
    S3 = S1 .* S2
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [0.0, 0.5, 0.0]
    
    S3 = @. S2 * (S2 + S1)
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [0.25, 0.75, 0.25]
    
    # SampledSpectrum with numbers
    
    S3 = S1 .* 2.0
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [0.0, 2.0, 0.0]
    
    S3 = S1 ./ 2.0
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [0.0, 0.5, 0.0]
    
    S3 = S1 .+ 2.0
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [2.0, 3.0, 2.0]
    
    S3 = 2.0 .* S1
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [0.0, 2.0, 0.0]
    
    
    S3 = 2.0 .+ S1
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [2.0, 3.0, 2.0]
    
    
    
    # SampledSpectrum with NoLight
    
    N = NoLight()
        
    S3 = N .* (S1 .+ S2)
    @test typeof(S3) == NoLight
    @test S3 == nolight
    
    S3 = N .+ (S1 .+ S2)
    @test typeof(S3) == SampledSpectrum
    @test S3.values ≈ [0.5, 1.5, 0.5]
    
end




end # main testset

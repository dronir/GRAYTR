
@testset "Samplers" begin

@testset "CameraSample" begin
    s = GRAYTR.CameraSample(4.0, 0.5, 0.75, 1.0)
    t = normalize(s, 2.0, 0.5)
    @test t.imgX ≈ 2.00
    @test t.imgY ≈ 1.0
    @test t.lensU ≈ 0.75
    @test t.lensV ≈ 1.0
end

@testset "StratifiedSampler" begin
    @testset "Creating" begin
        S = GRAYTR.StratifiedSampler(64, 64, 2)
        @test S.xstart == 1
        @test S.xend == 64
        @test S.ystart == 1
        @test S.yend == 64
        @test S.xs == 2
        @test S.ys == 2
        @test S.jitter == true
        @test S.xnorm == 64
        @test S.ynorm == 64
    end

    @testset "Subwindow computation" begin
        S = GRAYTR.StratifiedSampler(32, 64, 2)
        @test GRAYTR.compute_subwindow(S, 1, 4) == (1, 16, 1, 32)
        @test GRAYTR.compute_subwindow(S, 2, 4) == (17, 32, 1, 32)
        @test GRAYTR.compute_subwindow(S, 3, 4) == (1, 16, 33, 64)
        @test GRAYTR.compute_subwindow(S, 4, 4) == (17, 32, 33, 64)
        
        S = GRAYTR.StratifiedSampler(2, 2, 2)
        @test GRAYTR.compute_subwindow(S, 1, 4) == (1, 1, 1, 1)
        @test GRAYTR.compute_subwindow(S, 2, 4) == (2, 2, 1, 1)
        @test GRAYTR.compute_subwindow(S, 3, 4) == (1, 1, 2, 2)
        @test GRAYTR.compute_subwindow(S, 4, 4) == (2, 2, 2, 2)
        
        #S = GRAYTR.StratifiedSampler(2, 4, 2)
        #@test GRAYTR.compute_subwindow(S, 1, 4) == (1, 1, 1, 2)
        #@test GRAYTR.compute_subwindow(S, 2, 4) == (2, 2, 1, 2)
        #@test GRAYTR.compute_subwindow(S, 3, 4) == (1, 1, 3, 4)
        #@test GRAYTR.compute_subwindow(S, 4, 4) == (2, 2, 3, 4)
        #
        #S = GRAYTR.StratifiedSampler(4, 2, 2)
        #@test GRAYTR.compute_subwindow(S, 1, 4) == (1, 2, 1, 1)
        #@test GRAYTR.compute_subwindow(S, 2, 4) == (3, 4, 1, 1)
        #@test GRAYTR.compute_subwindow(S, 3, 4) == (1, 2, 2, 2)
        #@test GRAYTR.compute_subwindow(S, 4, 4) == (3, 4, 2, 2)
        
    end
    
    @testset "Subsampler" begin
        S = GRAYTR.StratifiedSampler(32, 64, 2)
        Sub = GRAYTR.get_subsampler(S, 2, 4)
        @test Sub.xstart == 17
        @test Sub.xend == 32
        @test Sub.ystart == 1
        @test Sub.yend == 32
        @test Sub.xs == 2
        @test Sub.ys == 2
        @test Sub.jitter == true
        @test Sub.xnorm == 32
        @test Sub.ynorm == 64
        
    end
    
    @testset "Finished state" begin
        S = GRAYTR.StratifiedSampler(4, 4, 2)
        @test !GRAYTR.finished(S, 2)
        @test !GRAYTR.finished(S, 15)
        @test GRAYTR.finished(S, 16)
        @test GRAYTR.finished(S, 17)
    end
    
    @testset "Getting samples" begin
        S = GRAYTR.StratifiedSampler(4, 4, 2)
        state = 0
        samples, state = GRAYTR.get_samples(S, state)
        
        @test length(samples) == 4
        @test eltype(samples) == GRAYTR.CameraSample
        @test state == 1
    end


end




end
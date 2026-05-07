using Test
using Aqua
using Documenter
using ExplicitImports
import RegisterMismatchCommon
using RegisterMismatchCommon
using RegisterCore: MismatchArray, NumDenom, maxshift, separate

# Tests are in RegisterMismatch and RegisterMismatchCuda

DocMeta.setdocmeta!(RegisterMismatchCommon, :DocTestSetup, :(using RegisterMismatchCommon); recursive=true)

@testset "RegisterMismatchCommon" begin
    @testset "Doctests" begin
        doctest(RegisterMismatchCommon; manual=false)
    end
    Aqua.test_all(RegisterMismatchCommon)
    @testset "ExplicitImports" begin
        test_explicit_imports(RegisterMismatchCommon)
    end

    @testset "set_FFTPROD" begin
        old = copy(RegisterMismatchCommon.FFTPROD)
        set_FFTPROD([2, 3, 5])
        @test RegisterMismatchCommon.FFTPROD == [2, 3, 5]
        set_FFTPROD(old)
        @test RegisterMismatchCommon.FFTPROD == old
    end

    @testset "tovec" begin
        v = [1, 2, 3]
        @test tovec(v) === v
        @test tovec((1, 2, 3)) == [1, 2, 3]
    end

    @testset "shiftrange" begin
        @test shiftrange(1:5, 3) == 4:8
        @test shiftrange(-2:2, -1) == -3:1
    end

    @testset "nanpad" begin
        a = [1.0 2.0; 3.0 4.0]
        b = [5.0 6.0; 7.0 8.0]
        # Same size: returns originals unchanged
        ap, bp = nanpad(a, b)
        @test ap === a
        @test bp === b

        # Different sizes: pads smaller to larger
        c = [1.0 2.0 3.0; 4.0 5.0 6.0]  # 2×3
        d = [7.0 8.0; 9.0 10.0]          # 2×2
        cp, dp = nanpad(c, d)
        @test size(cp) == size(dp) == (2, 3)
        @test cp[1, 1] == 1.0 && isnan(dp[1, 3])

        # Dimension mismatch
        @test_throws ErrorException nanpad([1.0], [1.0 2.0])
    end

    @testset "mismatch0 (fixed/moving)" begin
        a = [1.0 2.0; 3.0 4.0]
        # Identical arrays → num==0
        nd = mismatch0(a, a)
        @test nd.num == 0.0
        @test nd.denom > 0.0

        # Different arrays → positive num
        b = [2.0 3.0; 4.0 5.0]
        nd2 = mismatch0(a, b)
        @test nd2.num > 0.0

        # :pixels normalization
        nd3 = mismatch0(a, b; normalization = :pixels)
        @test nd3.denom == 4.0   # 4 finite pixels

        # NaN elements are skipped
        c = [1.0 NaN; 3.0 4.0]
        nd4 = mismatch0(a, c)
        @test nd4.denom ≈ a[1,1]^2 + c[1,1]^2 + a[2,1]^2 + c[2,1]^2 + a[2,2]^2 + c[2,2]^2

        # Size mismatch error
        @test_throws DimensionMismatch mismatch0(a, [1.0 2.0])

        # Unrecognized normalization
        @test_throws ErrorException mismatch0(a, a; normalization = :bogus)
    end

    @testset "mismatch0 (array-of-MismatchArrays)" begin
        mm = MismatchArray(Float64, 3, 3)
        mm[0, 0] = NumDenom{Float64}(1.0, 2.0)
        mms = [mm]
        nd = mismatch0(mms)
        @test nd == NumDenom{Float64}(1.0, 2.0)
    end

    @testset "aperture_grid" begin
        # 1D, single cell: center at midpoint
        g = aperture_grid((10,), (1,))
        @test size(g) == (1,)
        @test g[1] == (5.5,)

        # 1D, two cells: at corners
        g2 = aperture_grid((10,), (2,))
        @test g2[1] == (1.0,)
        @test g2[2] == (10.0,)

        # 2D grid
        g3 = aperture_grid((10, 20), (2, 3))
        @test size(g3) == (2, 3)
        @test g3[1, 1] == (1.0, 1.0)
        @test g3[2, 3] == (10.0, 20.0)

        # Length mismatch error
        @test_throws ErrorException aperture_grid((10, 20), (2,))
    end

    @testset "aperture_range" begin
        rng = aperture_range((5.0, 10.0), (4.0, 6.0))
        @test length(rng) == 2
        @test rng[1] == 3:6
        @test rng[2] == 7:12

        @test_throws ErrorException aperture_range((1.0,), (1.0, 2.0))
    end

    @testset "each_point" begin
        # Array-of-tuples
        pts = [(1.0, 2.0), (3.0, 4.0)]
        collected = collect(each_point(pts))
        @test collected == pts

        # Real matrix (columns are points)
        M = [1.0 3.0; 2.0 4.0]  # 2×2: two 2D points
        collected2 = collect(each_point(M))
        @test collected2[1] == [1.0, 2.0]
        @test collected2[2] == [3.0, 4.0]
    end

    @testset "checksize_maxshift" begin
        mm = MismatchArray(Float64, 3, 5)
        @test checksize_maxshift(mm, (1, 2)) === nothing

        # Wrong ndims
        @test_throws ErrorException checksize_maxshift(mm, (1,))

        # Wrong size along dim 1
        @test_throws ErrorException checksize_maxshift(mm, (2, 2))
    end

    @testset "padsize" begin
        # dim==1 uses nextpow(2,...)
        @test padsize(10, 3, 1) == nextpow(2, 16)
        # dim==2 uses nextprod
        @test padsize(10, 3, 2) == nextprod(RegisterMismatchCommon.FFTPROD, 16)
        # maxshift==0: no padding, no FFT
        @test padsize(10, 0, 1) == 10
        @test padsize(10, 0, 2) == 10

        # Dims version
        ps = padsize((8, 12), (2, 3))
        @test ps[1] == padsize(8, 2, 1)
        @test ps[2] == padsize(12, 3, 2)
    end

    @testset "padranges" begin
        rngs = padranges((10, 12), (2, 3))
        @test rngs[1].start == 1 - 2
        @test rngs[2].start == 1 - 3
        for (i, r) in enumerate(rngs)
            @test length(r) >= 10 + 2 * [2, 3][i]
        end
    end

    @testset "issamesize / assertsamesize" begin
        a = zeros(3, 4)
        b = zeros(3, 4)
        c = zeros(3, 5)

        @test RegisterMismatchCommon.issamesize(a, b)
        @test !RegisterMismatchCommon.issamesize(a, c)
        @test !RegisterMismatchCommon.issamesize(a, zeros(2, 3, 4))

        # Array vs. indices form
        @test RegisterMismatchCommon.issamesize(a, (1:3, 1:4))
        @test !RegisterMismatchCommon.issamesize(a, (1:3, 1:5))
        @test !RegisterMismatchCommon.issamesize(a, (1:3,))

        @test assertsamesize(a, b) === nothing
        @test_throws ErrorException assertsamesize(a, c)
    end

    @testset "allocate_mmarrays" begin
        # Array-of-tuples centers
        centers = [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]
        mms = allocate_mmarrays(Float32, centers, (2, 3))
        @test size(mms) == (3,)
        @test maxshift(mms[1]) == (2, 3)

        # Real matrix centers: rows are coordinates, columns are points
        # A 2×3 matrix means 3 points in 2D space → N=1 (ndims-1=1), 3 mmarrays
        M = [1.0 3.0 5.0; 2.0 4.0 6.0]
        mms2 = allocate_mmarrays(Float32, M, (2,))
        @test size(mms2) == (3,)

        # Gridsize tuple
        mms3 = allocate_mmarrays(Float32, (2, 3), (1, 2))
        @test size(mms3) == (2, 3)
        @test maxshift(mms3[1, 1]) == (1, 2)
    end

    @testset "truncatenoise!" begin
        mm = MismatchArray(Float64, 3, 3)
        for I in eachindex(mm)
            mm[I] = NumDenom{Float64}(1.0, Float64(abs(I[1]) + abs(I[2]) + 1))
        end
        thresh = 2.5
        truncatenoise!(mm, thresh)
        for I in eachindex(mm)
            d = abs(I[1]) + abs(I[2]) + 1
            if d <= thresh
                @test mm[I] == NumDenom{Float64}(0.0, 0.0)
            end
        end

        # Array-of-MismatchArrays form
        mms = [MismatchArray(Float64, 3, 3) for _ in 1:2]
        for mm in mms, I in eachindex(mm)
            mm[I] = NumDenom{Float64}(1.0, 0.5)
        end
        result = truncatenoise!(mms, 1.0)
        @test result === mms
        for mm in mms, I in eachindex(mm)
            @test mm[I] == NumDenom{Float64}(0.0, 0.0)
        end
    end

    @testset "correctbias!" begin
        # 2D MismatchArray: row 0 and column 0 are suspect
        mm = MismatchArray(Float64, 5, 5)  # shifts -2:2 in each dim
        good = NumDenom{Float64}(1.0, 1.0)
        bad  = NumDenom{Float64}(99.0, 99.0)
        for I in eachindex(mm)
            mm[I] = (I[1] == 0 || I[2] == 0) ? bad : good
        end
        correctbias!(mm)
        # All suspect entries should now be replaced
        for I in eachindex(mm)
            if I[1] == 0 || I[2] == 0
                @test mm[I] != bad
            end
        end

        # Array-of-MismatchArrays form
        mms = [MismatchArray(Float64, 5, 5) for _ in 1:2]
        for mm2 in mms, I in eachindex(mm2)
            mm2[I] = (I[1] == 0 || I[2] == 0) ? bad : good
        end
        result = correctbias!(mms)
        @test result === mms
        for mm2 in mms, I in eachindex(mm2)
            if I[1] == 0 || I[2] == 0
                @test mm2[I] != bad
            end
        end
    end

    @testset "default_aperture_width" begin
        img = zeros(Float32, 20, 30)
        # Single cell: full image width
        w = default_aperture_width(img, (1, 1))
        @test w == (20.0, 30.0)

        # 2×3 grid
        w2 = default_aperture_width(img, (2, 3))
        @test length(w2) == 2

        # gridsize too large
        @test_throws ErrorException default_aperture_width(img, (25, 3))
    end

    @testset "nanval" begin
        @test isnan(RegisterMismatchCommon.nanval(Float64))
        @test isnan(RegisterMismatchCommon.nanval(Float32))
        # Non-float falls back to Float32 NaN
        v = RegisterMismatchCommon.nanval(Int32)
        @test isa(v, Float32) && isnan(v)
    end

    @testset "computeoverlap" begin
        # 2×2 grid on a 10×10 image with 6×6 blocksize → spacing = ceil(10/1) = 10, overlap = 6-10 = -4
        ov = RegisterMismatchCommon.computeoverlap((10, 10), (6, 6), (2, 2))
        @test length(ov) == 2
        @test ov == [6 - 10, 6 - 10]
        # gridsize 1 along first dim: gsz1 = max.(1, [0, 1]) = [1, 1]
        ov2 = RegisterMismatchCommon.computeoverlap((10, 10), (6, 6), (1, 2))
        @test length(ov2) == 2
    end
end

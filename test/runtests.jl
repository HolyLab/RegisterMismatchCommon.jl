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

    @testset "FFTPROD" begin
        old = copy(RegisterMismatchCommon.FFTPROD[])
        RegisterMismatchCommon.FFTPROD[] = [2, 3, 5]
        @test RegisterMismatchCommon.FFTPROD[] == [2, 3, 5]
        RegisterMismatchCommon.FFTPROD[] = old
        @test RegisterMismatchCommon.FFTPROD[] == old
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

    @testset "mismatch_zeroshift (fixed/moving)" begin
        a = [1.0 2.0; 3.0 4.0]
        # Identical arrays → num==0
        nd = mismatch_zeroshift(a, a)
        @test nd.num == 0.0
        @test nd.denom > 0.0

        # Different arrays → positive num
        b = [2.0 3.0; 4.0 5.0]
        nd2 = mismatch_zeroshift(a, b)
        @test nd2.num > 0.0

        # :pixels normalization
        nd3 = mismatch_zeroshift(a, b; normalization = :pixels)
        @test nd3.denom == 4.0   # 4 finite pixels

        # NaN elements are skipped
        c = [1.0 NaN; 3.0 4.0]
        nd4 = mismatch_zeroshift(a, c)
        @test nd4.denom ≈ a[1,1]^2 + c[1,1]^2 + a[2,1]^2 + c[2,1]^2 + a[2,2]^2 + c[2,2]^2

        # Size mismatch error
        @test_throws DimensionMismatch mismatch_zeroshift(a, [1.0 2.0])

        # Unrecognized normalization
        @test_throws ErrorException mismatch_zeroshift(a, a; normalization = :bogus)

        # Deprecated alias still works
        @test_deprecated mismatch0(a, a)
    end

    @testset "mismatch_zeroshift (array-of-MismatchArrays)" begin
        mm = MismatchArray(Float64, 3, 3)
        mm[0, 0] = NumDenom{Float64}(1.0, 2.0)
        mms = [mm]
        nd = mismatch_zeroshift(mms)
        @test nd == NumDenom{Float64}(1.0, 2.0)

        # Deprecated alias still works
        @test_deprecated mismatch0(mms)
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
        @test padsize(10, 3, 2) == nextprod(RegisterMismatchCommon.FFTPROD[], 16)
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

    @testset "issamesize / checksamesize" begin
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

        @test checksamesize(a, b) === nothing
        @test_throws ErrorException checksamesize(a, c)

        # Deprecated alias still works (forwards to checksamesize)
        @test_deprecated assertsamesize(a, b)
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

    @testset "truncatenoise (non-mutating)" begin
        mm = MismatchArray(Float64, 3, 3)
        for I in eachindex(mm)
            mm[I] = NumDenom{Float64}(1.0, Float64(abs(I[1]) + abs(I[2]) + 1))
        end
        thresh = 2.5
        mm_orig = copy(mm)
        mm_trunc = truncatenoise(mm, thresh)
        @test mm_trunc !== mm
        # Original untouched
        for I in eachindex(mm)
            @test mm[I] == mm_orig[I]
        end
        # Result matches in-place version
        mm_inplace = truncatenoise!(copy(mm), thresh)
        for I in eachindex(mm)
            @test mm_trunc[I] == mm_inplace[I]
        end

        # Array-of-MismatchArrays form
        mms = [MismatchArray(Float64, 3, 3) for _ in 1:2]
        for mm2 in mms, I in eachindex(mm2)
            mm2[I] = NumDenom{Float64}(1.0, 0.5)
        end
        mms_orig = deepcopy(mms)
        mms_trunc = truncatenoise(mms, 1.0)
        @test mms_trunc !== mms
        # Originals untouched
        for k in eachindex(mms), I in eachindex(mms[k])
            @test mms[k][I] == mms_orig[k][I]
        end
        # Result was thresholded
        for mm2 in mms_trunc, I in eachindex(mm2)
            @test mm2[I] == NumDenom{Float64}(0.0, 0.0)
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

    @testset "correctbias (non-mutating)" begin
        mm = MismatchArray(Float64, 5, 5)
        good = NumDenom{Float64}(1.0, 1.0)
        bad  = NumDenom{Float64}(99.0, 99.0)
        for I in eachindex(mm)
            mm[I] = (I[1] == 0 || I[2] == 0) ? bad : good
        end
        mm_orig = copy(mm)
        mm_corr = correctbias(mm)
        @test mm_corr !== mm
        # Original untouched
        for I in eachindex(mm)
            @test mm[I] == mm_orig[I]
        end
        # Result matches in-place version
        mm_inplace = correctbias!(copy(mm))
        for I in eachindex(mm)
            @test mm_corr[I] == mm_inplace[I]
        end

        # Array-of-MismatchArrays form
        mms = [MismatchArray(Float64, 5, 5) for _ in 1:2]
        for mm2 in mms, I in eachindex(mm2)
            mm2[I] = (I[1] == 0 || I[2] == 0) ? bad : good
        end
        mms_orig = deepcopy(mms)
        mms_corr = correctbias(mms)
        @test mms_corr !== mms
        # Originals untouched
        for k in eachindex(mms), I in eachindex(mms[k])
            @test mms[k][I] == mms_orig[k][I]
        end
        # Result matches in-place version
        for k in eachindex(mms_corr), I in eachindex(mms_corr[k])
            if I[1] == 0 || I[2] == 0
                @test mms_corr[k][I] != bad
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

    @testset "mismatch / mismatch_apertures (promoting dispatch)" begin
        # General protocol stubs — lower specificity than the Matrix{Float32} stub
        # defined in the register_translate testset below, so both coexist without conflict.
        function RegisterMismatchCommon.mismatch(::Type{T}, f, m, ms; kwargs...) where T
            MismatchArray(T, (2 .* collect(ms) .+ 1)...)
        end
        function RegisterMismatchCommon.mismatch_apertures(::Type{T}, f, m,
                                                            aperture_centers,
                                                            aperture_width, ms;
                                                            kwargs...) where T
            allocate_mmarrays(T, aperture_centers, ms)
        end

        # T1: non-AbstractFloat arrays promote to Float32 (line 80)
        result = mismatch(ones(Int, 4, 4), ones(Int, 4, 4), (1, 1))
        @test eltype(result) === NumDenom{Float32}

        # T2: mismatch_apertures promoting dispatch
        centers = [(4.0, 4.0), (8.0, 4.0)]
        width   = (4.0, 4.0)
        # AbstractFloat array → dispatches typed (line 82); Float64 stays Float64
        mms64 = mismatch_apertures(ones(Float64, 12, 8), ones(Float64, 12, 8),
                                   centers, width, (1, 1))
        @test eltype(first(mms64)) === NumDenom{Float64}
        # Non-Float array → promoted to Float32 (line 83)
        mms_int = mismatch_apertures(ones(Int, 12, 8), ones(Int, 12, 8),
                                     centers, width, (1, 1))
        @test eltype(first(mms_int)) === NumDenom{Float32}

        # T3: mismatch_apertures gridsize form — builds aperture_centers internally (lines 85-89)
        fixed = ones(Float32, 16, 16)
        mms_gs = mismatch_apertures(Float32, fixed, fixed, (2, 2), (1, 1))
        @test size(mms_gs) == (2, 2)
        @test eltype(first(mms_gs)) === NumDenom{Float32}
    end

    @testset "register_translate (thresh shim)" begin
        # Stub out the protocol mismatch method so register_translate can be exercised
        # without a concrete RegisterMismatch package.
        # MismatchArray(T, n, m) takes total size; for maxshift (k,k), size is (2k+1, 2k+1).
        function RegisterMismatchCommon.mismatch(::Type{Float32},
                                                  fixed::Matrix{Float32},
                                                  moving::Matrix{Float32},
                                                  maxshift::Dims{2}; kwargs...)
            sz = 2 .* maxshift .+ 1
            mm = MismatchArray(Float32, sz...)
            mm[0, 0] = NumDenom{Float32}(0.0f0, 10.0f0)   # clear winner: zero mismatch, high denom
            mm[1, 0] = NumDenom{Float32}(5.0f0, 10.0f0)
            return mm
        end

        fixed  = zeros(Float32, 4, 4)
        moving = zeros(Float32, 4, 4)

        # keyword form returns the expected shift
        result_kw = register_translate(fixed, moving, (1, 1); thresh=0.0)
        @test result_kw == CartesianIndex(0, 0)

        # default thresh (no keyword) also works
        result_default = register_translate(fixed, moving, (1, 1))
        @test result_default == CartesianIndex(0, 0)

        # positional thresh form is no longer supported (removed in v1.0.2)
        @test_throws MethodError register_translate(fixed, moving, (1, 1), 0.0)
    end
end

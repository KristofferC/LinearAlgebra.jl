module TestUnitfulLinAlg

using Test, LinearAlgebra, Random

Random.seed!(1234321)

const BASE_TEST_PATH = joinpath(Sys.BINDIR, "..", "share", "julia", "test")
isdefined(Main, :Furlongs) || @eval Main include(joinpath($(BASE_TEST_PATH), "testhelpers", "Furlongs.jl"))
using .Main.Furlongs

LinearAlgebra.sylvester(a::Furlong,b::Furlong,c::Furlong) = -c / (a + b)

getval(x) = x
getval(x::Furlong) = x.val

# specialized unitful multiplication/solves
@testset "unitful specialized mul/div" begin
    n = 10 #Size of test matrix
    function _bidiagdivmultest(T,
        x,
        typemul=T.uplo == 'U' ? UpperTriangular : Matrix,
        typediv=T.uplo == 'U' ? UpperTriangular : Matrix,
        typediv2=T.uplo == 'U' ? UpperTriangular : Matrix)
        TM = Matrix(T)
        @test map(getval, (T * x)::typemul) ≈ map(getval, TM * x)
        @test map(getval, (x * T)::typemul) ≈ map(getval, x * TM)
        @test map(getval, (x \ T)::typediv) ≈ map(getval, x \ TM)
        @test map(getval, (T / x)::typediv) ≈ map(getval, TM / x)
        if !isa(x, Number)
            @test map(getval, Array((T \ x)::typediv2)) ≈ map(getval, Array(TM \ x))
            @test map(getval, Array((x / T)::typediv2)) ≈ map(getval, Array(x / TM))
        end
        return nothing
    end
    @testset for relty in (Int, Float64, BigFloat), elty in (relty, Complex{relty})
        if relty <: AbstractFloat
            dv = convert(Vector{elty}, randn(n))
            ev = convert(Vector{elty}, randn(n - 1))
            if (elty <: Complex)
                dv += im * convert(Vector{elty}, randn(n))
                ev += im * convert(Vector{elty}, randn(n - 1))
            end
        elseif relty <: Integer
            dv = convert(Vector{elty}, rand(1:10, n))
            ev = convert(Vector{elty}, rand(1:10, n - 1))
            if (elty <: Complex)
                dv += im * convert(Vector{elty}, rand(1:10, n))
                ev += im * convert(Vector{elty}, rand(1:10, n - 1))
            end
        end
        @testset for uplo in (:U, :L)
            T = Bidiagonal(dv, ev, uplo)
            A = Matrix(T)
            for t in (T, Furlong.(T)), (A, dv, ev) in ((A, dv, ev), (Furlong.(A), Furlong.(dv), Furlong.(ev)))
                any(x -> x <: Furlong, (eltype(t), eltype(A))) || continue
                _bidiagdivmultest(t, 5, Bidiagonal, Bidiagonal)
                _bidiagdivmultest(t, 5I, Bidiagonal, Bidiagonal, t.uplo == 'U' ? UpperTriangular : LowerTriangular)
                _bidiagdivmultest(t, Diagonal(dv), Bidiagonal, Bidiagonal, t.uplo == 'U' ? UpperTriangular : LowerTriangular)
                _bidiagdivmultest(t, UpperTriangular(A))
                _bidiagdivmultest(t, UnitUpperTriangular(A))
                _bidiagdivmultest(t, LowerTriangular(A), t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix)
                _bidiagdivmultest(t, UnitLowerTriangular(A), t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix, t.uplo == 'L' ? LowerTriangular : Matrix)
                _bidiagdivmultest(t, Bidiagonal(dv, ev, :U), Matrix, Matrix, Matrix)
                _bidiagdivmultest(t, Bidiagonal(dv, ev, :L), Matrix, Matrix, Matrix)
            end
        end
    end
end

# diagonal
@testset for relty in (Float32, Float64, BigFloat), elty in (relty, Complex{relty})
    n = 12
    dd = convert(Vector{elty}, randn(n))
    D = Diagonal(dd)
    @testset "svd/eigen with Diagonal{Furlong}" begin
        Du = Furlong.(D)
        @test Du isa Diagonal{<:Furlong{1}}
        F = svd(Du)
        U, s, V = F
        @test map(getval, Matrix(F)) ≈ map(getval, Du)
        @test svdvals(Du) == s
        @test U isa AbstractMatrix{<:Furlong{0}}
        @test V isa AbstractMatrix{<:Furlong{0}}
        @test s isa AbstractVector{<:Furlong{1}}
        E = eigen(Du)
        vals, vecs = E
        @test Matrix(E) == Du
        @test vals isa AbstractVector{<:Furlong{1}}
        @test vecs isa AbstractMatrix{<:Furlong{0}}
    end
end

# givens
@testset "testing dimensions with Furlongs #36430" begin
    @test_throws MethodError givens(Furlong(1.0), Furlong(2.0), 1, 2)
end

# hessenberg
@testset "dimensional Hessenberg" begin
    n = 10
    Random.seed!(1234321)

    Areal = randn(n, n) / 2
    @testset "Preserve UpperHessenberg shape (issue #39388)" begin
        H = UpperHessenberg(Furlong.(Areal))
        A = Furlong.(rand(n, n))
        d = Furlong.(rand(n))
        dl = Furlong.(rand(n - 1))
        du = Furlong.(rand(n - 1))
        us = Furlong(1) * I

        @testset "$op" for op = (+, -)
            for x = (us, Diagonal(d), Bidiagonal(d, dl, :U), Bidiagonal(d, dl, :L),
                Tridiagonal(dl, d, du), SymTridiagonal(d, dl),
                UpperTriangular(A), UnitUpperTriangular(A))
                @test op(H, x) == op(Array(H), x)
                @test op(x, H) == op(x, Array(H))
                @test op(H, x) isa UpperHessenberg
                @test op(x, H) isa UpperHessenberg
            end
        end
    end
    H = UpperHessenberg(Furlong.(Areal))
    A = randn(n, n)
    d = randn(n)
    dl = randn(n - 1)
    for A in (A, Furlong.(A))
        @testset "Multiplication/division Furlong" begin
            for x = (5, 5I, Diagonal(d), Bidiagonal(d, dl, :U),
                UpperTriangular(A), UnitUpperTriangular(A))
                @test map(getval, (H * x)::UpperHessenberg) ≈ map(getval, Array(H) * x)
                @test map(getval, (x * H)::UpperHessenberg) ≈ map(getval, x * Array(H))
                @test map(getval, (H / x)::UpperHessenberg) ≈ map(getval, Array(H) / x)
                @test map(getval, (x \ H)::UpperHessenberg) ≈ map(getval, x \ Array(H))
            end
            x = Bidiagonal(d, dl, :L)
            @test H * x == Array(H) * x
            @test x * H == x * Array(H)
            @test H / x == Array(H) / x
            @test x \ H == x \ Array(H)
        end
    end
end

# lu
@testset "lu factorization with dimension type" begin
    n = 4
    A = Matrix(Furlong(1.0) * I, n, n)
    F = lu(A).factors
    @test Diagonal(F) == Diagonal(A)
    # upper triangular part has a unit Furlong{1}
    @test all(x -> typeof(x) == Furlong{1,Float64}, F[i, j] for j = 1:n for i = 1:j)
    # lower triangular part is unitless Furlong{0}
    @test all(x -> typeof(x) == Furlong{0,Float64}, F[i, j] for j = 1:n for i = j+1:n)
end

# special
@testset "zero and one for unitful structured matrices" begin
    # eltype with dimensions
    D0 = Diagonal{Furlong{0,Int64}}([1, 2, 3, 4])
    Bu0 = Bidiagonal{Furlong{0,Int64}}([1, 2, 3, 4], [1, 2, 3], 'U')
    Bl0 = Bidiagonal{Furlong{0,Int64}}([1, 2, 3, 4], [1, 2, 3], 'L')
    T0 = Tridiagonal{Furlong{0,Int64}}([1, 2, 3], [1, 2, 3, 4], [1, 2, 3])
    S0 = SymTridiagonal{Furlong{0,Int64}}([1, 2, 3, 4], [1, 2, 3])
    F2 = Furlongs.Furlong{2}(1)
    D2 = Diagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2)
    Bu2 = Bidiagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2, [1, 2, 3] .* F2, 'U')
    Bl2 = Bidiagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2, [1, 2, 3] .* F2, 'L')
    T2 = Tridiagonal{Furlong{2,Int64}}([1, 2, 3] .* F2, [1, 2, 3, 4] .* F2, [1, 2, 3] .* F2)
    S2 = SymTridiagonal{Furlong{2,Int64}}([1, 2, 3, 4] .* F2, [1, 2, 3] .* F2)
    mats = Any[D0, Bu0, Bl0, T0, S0, D2, Bu2, Bl2, T2, S2]
    for A in mats
        @test iszero(zero(A))
        @test isone(one(A))
        @test zero(A) == zero(Matrix(A))
        @test one(A) == one(Matrix(A))
        @test eltype(one(A)) == typeof(one(eltype(A)))
    end
end

# triangular
@testset "triangular: dimensional correctness" begin
    A = UpperTriangular([Furlong(1) Furlong(4); Furlong(0) Furlong(1)])
    @test sqrt(A)::UpperTriangular == Furlong{1 // 2}.(UpperTriangular([1 2; 0 1]))
    @test inv(A)::UpperTriangular == Furlong{-1}.(UpperTriangular([1 -4; 0 1]))
    B = UnitUpperTriangular([Furlong(1) Furlong(4); Furlong(0) Furlong(1)])
    @test sqrt(B)::UnitUpperTriangular == Furlong{1 // 2}.(UpperTriangular([1 2; 0 1]))
    @test inv(B)::UnitUpperTriangular == Furlong{-1}.(UpperTriangular([1 -4; 0 1]))
    b = [Furlong(5), Furlong(8)]
    @test (A \ b)::Vector{<:Furlong{0}} == (B \ b)::Vector{<:Furlong{0}} == Furlong{0}.([-27, 8])
    C = LowerTriangular([Furlong(1) Furlong(0); Furlong(4) Furlong(1)])
    @test sqrt(C)::LowerTriangular == Furlong{1 // 2}.(LowerTriangular([1 0; 2 1]))
    @test inv(C)::LowerTriangular == Furlong{-1}.(LowerTriangular([1 0; -4 1]))
    D = UnitLowerTriangular([Furlong(1) Furlong(0); Furlong(4) Furlong(1)])
    @test sqrt(D)::UnitLowerTriangular == Furlong{1 // 2}.(UnitLowerTriangular([1 0; 2 1]))
    @test inv(D)::UnitLowerTriangular == Furlong{-1}.(UnitLowerTriangular([1 0; -4 1]))
    b = [Furlong(5), Furlong(8)]
    @test (C \ b)::Vector{<:Furlong{0}} == (D \ b)::Vector{<:Furlong{0}} == Furlong{0}.([5, -12])
end

end # module TestUnitfulLinAlg

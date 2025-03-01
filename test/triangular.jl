# This file is a part of Julia. License is MIT: https://julialang.org/license

module TestTriangular

using Test, LinearAlgebra, Random
using LinearAlgebra: errorbounds, transpose!, BandIndex

const BASE_TEST_PATH = joinpath(Sys.BINDIR, "..", "share", "julia", "test")

isdefined(Main, :SizedArrays) || @eval Main include(joinpath($(BASE_TEST_PATH), "testhelpers", "SizedArrays.jl"))
using .Main.SizedArrays

isdefined(Main, :FillArrays) || @eval Main include(joinpath($(BASE_TEST_PATH), "testhelpers", "FillArrays.jl"))
using .Main.FillArrays

n = 9
Random.seed!(123)

struct MyTriangular{T, A<:LinearAlgebra.AbstractTriangular{T}} <: LinearAlgebra.AbstractTriangular{T}
    data :: A
end
Base.size(A::MyTriangular) = size(A.data)
Base.getindex(A::MyTriangular, i::Int, j::Int) = A.data[i,j]

@testset "Test basic type functionality" begin
    @test_throws DimensionMismatch LowerTriangular(randn(5, 4))
    @test LowerTriangular(randn(3, 3)) |> t -> [size(t, i) for i = 1:3] == [size(Matrix(t), i) for i = 1:3]
end

@testset "non-strided arithmetic" begin
    for (T,T1) in ((UpperTriangular, UnitUpperTriangular), (LowerTriangular, UnitLowerTriangular))
        U = T(reshape(1:16, 4, 4))
        M = Matrix(U)
        @test -U == -M
        U1 = T1(reshape(1:16, 4, 4))
        M1 = Matrix(U1)
        @test -U1 == -M1
        for op in (+, -)
            for (A, MA) in ((U, M), (U1, M1)), (B, MB) in ((U, M), (U1, M1))
                @test op(A, B) == op(MA, MB)
            end
        end
        @test imag(U) == zero(U)
    end
end

@testset "Matrix square root" begin
    Atn = UpperTriangular([-1 1 2; 0 -2 2; 0 0 -3])
    Atp = UpperTriangular([1 1 2; 0 2 2; 0 0 3])
    Atu = UnitUpperTriangular([1 1 2; 0 1 2; 0 0 1])
    @test sqrt(Atn) |> t->t*t ≈ Atn
    @test sqrt(Atn) isa UpperTriangular
    @test typeof(sqrt(Atn)[1,1]) <: Complex
    @test sqrt(Atp) |> t->t*t ≈ Atp
    @test sqrt(Atp) isa UpperTriangular
    @test typeof(sqrt(Atp)[1,1]) <: Real
    @test typeof(sqrt(complex(Atp))[1,1]) <: Complex
    @test sqrt(Atu) |> t->t*t ≈ Atu
    @test sqrt(Atu) isa UnitUpperTriangular
    @test typeof(sqrt(Atu)[1,1]) <: Real
    @test typeof(sqrt(complex(Atu))[1,1]) <: Complex
end

@testset "matrix square root quasi-triangular blockwise" begin
    @testset for T in (Float32, Float64, ComplexF32, ComplexF64)
        A = schur(rand(T, 100, 100)^2).T
        @test LinearAlgebra.sqrt_quasitriu(A; blockwidth=16)^2 ≈ A
    end
    n = 256
    A = rand(ComplexF64, n, n)
    U = schur(A).T
    Ubig = Complex{BigFloat}.(U)
    @test LinearAlgebra.sqrt_quasitriu(U; blockwidth=64) ≈ LinearAlgebra.sqrt_quasitriu(Ubig; blockwidth=64)
end

@testset "sylvester quasi-triangular blockwise" begin
    @testset for T in (Float32, Float64, ComplexF32, ComplexF64), m in (15, 40), n in (15, 45)
        A = schur(rand(T, m, m)).T
        B = schur(rand(T, n, n)).T
        C = randn(T, m, n)
        Ccopy = copy(C)
        X = LinearAlgebra._sylvester_quasitriu!(A, B, C; blockwidth=16)
        @test X === C
        @test A * X + X * B ≈ -Ccopy

        @testset "test raise=false does not break recursion" begin
            Az = zero(A)
            Bz = zero(B)
            C2 = copy(Ccopy)
            @test_throws LAPACKException LinearAlgebra._sylvester_quasitriu!(Az, Bz, C2; blockwidth=16)
            m == n || @test any(C2 .== Ccopy)  # recursion broken
            C3 = copy(Ccopy)
            X3 = LinearAlgebra._sylvester_quasitriu!(Az, Bz, C3; blockwidth=16, raise=false)
            @test !any(X3 .== Ccopy)  # recursion not broken
        end
    end
end

@testset "check matrix logarithm type-inferable" for elty in (Float32,Float64,ComplexF32,ComplexF64)
    A = UpperTriangular(exp(triu(randn(elty, n, n))))
    @inferred Union{typeof(A),typeof(complex(A))} log(A)
    @test exp(Matrix(log(A))) ≈ A
    if elty <: Real
        @test typeof(log(A)) <: UpperTriangular{elty}
        @test typeof(log(complex(A))) <: UpperTriangular{complex(elty)}
        @test isreal(log(complex(A)))
        @test log(complex(A)) ≈ log(A)
    end

    Au = UnitUpperTriangular(exp(triu(randn(elty, n, n), 1)))
    @inferred Union{typeof(A),typeof(complex(A))} log(Au)
    @test exp(Matrix(log(Au))) ≈ Au
    if elty <: Real
        @test typeof(log(Au)) <: UpperTriangular{elty}
        @test typeof(log(complex(Au))) <: UpperTriangular{complex(elty)}
        @test isreal(log(complex(Au)))
        @test log(complex(Au)) ≈ log(Au)
    end
end

Areal   = randn(n, n)/2
Aimg    = randn(n, n)/2
A2real  = randn(n, n)/2
A2img   = randn(n, n)/2

@testset for eltya in (Float32, Float64, ComplexF32, ComplexF64, BigFloat, Int)
    A = eltya == Int ? rand(1:7, n, n) : convert(Matrix{eltya}, eltya <: Complex ? complex.(Areal, Aimg) : Areal)
    εa = eps(abs(float(one(eltya))))

    for eltyb in (Float32, Float64, ComplexF32, ComplexF64)
        εb = eps(abs(float(one(eltyb))))
        ε = max(εa,εb)

        @testset "Solve upper triangular system" begin
            Atri = UpperTriangular(lu(A).U) |> t -> eltya <: Complex && eltyb <: Real ? real(t) : t # Here the triangular matrix can't be too badly conditioned
            b = convert(Matrix{eltyb}, Matrix(Atri)*fill(1., n, 2))
            x = Atri \ b

            # Test error estimates
            if eltya != BigFloat && eltyb != BigFloat
                for i = 1:2
                    @test norm(x[:,1] .- 1) <= errorbounds(UpperTriangular(A), x, b)[1][i]
                end
            end

            # Test forward error [JIN 5705] if this is not a BigFloat
            γ = n*ε/(1 - n*ε)
            if eltya != BigFloat
                bigA = big.(Atri)
                x̂ = fill(1., n, 2)
                for i = 1:size(b, 2)
                    @test norm(x̂[:,i] - x[:,i], Inf)/norm(x̂[:,i], Inf) <= condskeel(bigA, x̂[:,i])*γ/(1 - condskeel(bigA)*γ)
                end
            end

            # Test backward error [JIN 5705]
            for i = 1:size(b, 2)
                @test norm(abs.(b[:,i] - Atri*x[:,i]), Inf) <= γ * norm(Atri, Inf) * norm(x[:,i], Inf)
            end
        end

        @testset "Solve lower triangular system" begin
            Atri = LowerTriangular(lu(A).L) |> t -> eltya <: Complex && eltyb <: Real ? real(t) : t # Here the triangular matrix can't be too badly conditioned
            b = convert(Matrix{eltyb}, Matrix(Atri)*fill(1., n, 2))
            x = Atri \ b

            # Test error estimates
            if eltya != BigFloat && eltyb != BigFloat
                for i = 1:2
                    @test norm(x[:,1] .- 1) <= errorbounds(LowerTriangular(A), x, b)[1][i]
                end
            end

            # Test forward error [JIN 5705] if this is not a BigFloat
            γ = n*ε/(1 - n*ε)
            if eltya != BigFloat
                bigA = big.(Atri)
                x̂ = fill(1., n, 2)
                for i = 1:size(b, 2)
                    @test norm(x̂[:,i] - x[:,i], Inf)/norm(x̂[:,i], Inf) <= condskeel(bigA, x̂[:,i])*γ/(1 - condskeel(bigA)*γ)
                end
            end

            # Test backward error [JIN 5705]
            for i = 1:size(b, 2)
                @test norm(abs.(b[:,i] - Atri*x[:,i]), Inf) <= γ * norm(Atri, Inf) * norm(x[:,i], Inf)
            end
        end
    end
end

@testset "triangularity/diagonality of triangular views (#10742)" begin
    @test istril(UpperTriangular(diagm(0 => [1,2,3,4])))
    @test istriu(LowerTriangular(diagm(0 => [1,2,3,4])))
    @test isdiag(UpperTriangular(diagm(0 => [1,2,3,4])))
    @test isdiag(LowerTriangular(diagm(0 => [1,2,3,4])))
    @test !isdiag(UpperTriangular(rand(4, 4)))
    @test !isdiag(LowerTriangular(rand(4, 4)))
end

# Test throwing in fallbacks for non BlasFloat/BlasComplex in A_rdiv_Bx!
let n = 5
    A = rand(Float16, n, n)
    B = rand(Float16, n-1, n-1)
    @test_throws DimensionMismatch rdiv!(A, LowerTriangular(B))
    @test_throws DimensionMismatch rdiv!(A, UpperTriangular(B))
    @test_throws DimensionMismatch rdiv!(A, UnitLowerTriangular(B))
    @test_throws DimensionMismatch rdiv!(A, UnitUpperTriangular(B))

    @test_throws DimensionMismatch rdiv!(A, adjoint(LowerTriangular(B)))
    @test_throws DimensionMismatch rdiv!(A, adjoint(UpperTriangular(B)))
    @test_throws DimensionMismatch rdiv!(A, adjoint(UnitLowerTriangular(B)))
    @test_throws DimensionMismatch rdiv!(A, adjoint(UnitUpperTriangular(B)))

    @test_throws DimensionMismatch rdiv!(A, transpose(LowerTriangular(B)))
    @test_throws DimensionMismatch rdiv!(A, transpose(UpperTriangular(B)))
    @test_throws DimensionMismatch rdiv!(A, transpose(UnitLowerTriangular(B)))
    @test_throws DimensionMismatch rdiv!(A, transpose(UnitUpperTriangular(B)))
end

@test isdiag(LowerTriangular(UpperTriangular(randn(3,3))))
@test isdiag(UpperTriangular(LowerTriangular(randn(3,3))))

# Issue 16196
@test UpperTriangular(Matrix(1.0I, 3, 3)) \ view(fill(1., 3), [1,2,3]) == fill(1., 3)

@testset "reverse" begin
    A = randn(5, 5)
    for (T, Trev) in ((UpperTriangular, LowerTriangular),
            (UnitUpperTriangular, UnitLowerTriangular),
            (LowerTriangular, UpperTriangular),
            (UnitLowerTriangular, UnitUpperTriangular))
        A = T(randn(5, 5))
        AM = Matrix(A)
        @test reverse(A, dims=1) == reverse(AM, dims=1)
        @test reverse(A, dims=2) == reverse(AM, dims=2)
        @test reverse(A)::Trev == reverse(AM)
    end
end

# dimensional correctness:
const BASE_TEST_PATH = joinpath(Sys.BINDIR, "..", "share", "julia", "test")

isdefined(Main, :ImmutableArrays) || @eval Main include(joinpath($(BASE_TEST_PATH), "testhelpers", "ImmutableArrays.jl"))
using .Main.ImmutableArrays

@testset "AbstractArray constructor should preserve underlying storage type" begin
    # tests corresponding to #34995
    local m = 4
    local T, S = Float32, Float64
    immutablemat = ImmutableArray(randn(T,m,m))
    for TriType in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
        trimat = TriType(immutablemat)
        @test convert(AbstractArray{S}, trimat).data isa ImmutableArray{S}
        @test convert(AbstractMatrix{S}, trimat).data isa ImmutableArray{S}
        @test AbstractArray{S}(trimat).data isa ImmutableArray{S}
        @test AbstractMatrix{S}(trimat).data isa ImmutableArray{S}
        @test convert(AbstractArray{S}, trimat) == trimat
        @test convert(AbstractMatrix{S}, trimat) == trimat
    end
end

@testset "inplace mul of appropriate types should preserve triagular structure" begin
    for elty1 in (Float64, ComplexF32), elty2 in (Float64, ComplexF32)
        T = promote_type(elty1, elty2)
        M1 = rand(elty1, 5, 5)
        M2 = rand(elty2, 5, 5)
        A = UpperTriangular(M1)
        A2 = UpperTriangular(M2)
        Au = UnitUpperTriangular(M1)
        Au2 = UnitUpperTriangular(M2)
        B = LowerTriangular(M1)
        B2 = LowerTriangular(M2)
        Bu = UnitLowerTriangular(M1)
        Bu2 = UnitLowerTriangular(M2)

        @test mul!(similar(A), A, A)::typeof(A) == A*A
        @test mul!(similar(A, T), A, A2) ≈ A*A2
        @test mul!(similar(A, T), A2, A) ≈ A2*A
        @test mul!(typeof(similar(A, T))(A), A, A2, 2.0, 3.0) ≈ 2.0*A*A2 + 3.0*A
        @test mul!(typeof(similar(A2, T))(A2), A2, A, 2.0, 3.0) ≈ 2.0*A2*A + 3.0*A2

        @test mul!(similar(A), A, Au)::typeof(A) == A*Au
        @test mul!(similar(A), Au, A)::typeof(A) == Au*A
        @test mul!(similar(Au), Au, Au)::typeof(Au) == Au*Au
        @test mul!(similar(A, T), A, Au2) ≈ A*Au2
        @test mul!(similar(A, T), Au2, A) ≈ Au2*A
        @test mul!(similar(Au2), Au2, Au2) == Au2*Au2

        @test mul!(similar(B), B, B)::typeof(B) == B*B
        @test mul!(similar(B, T), B, B2) ≈ B*B2
        @test mul!(similar(B, T), B2, B) ≈ B2*B
        @test mul!(typeof(similar(B, T))(B), B, B2, 2.0, 3.0) ≈ 2.0*B*B2 + 3.0*B
        @test mul!(typeof(similar(B2, T))(B2), B2, B, 2.0, 3.0) ≈ 2.0*B2*B + 3.0*B2

        @test mul!(similar(B), B, Bu)::typeof(B) == B*Bu
        @test mul!(similar(B), Bu, B)::typeof(B) == Bu*B
        @test mul!(similar(Bu), Bu, Bu)::typeof(Bu) == Bu*Bu
        @test mul!(similar(B, T), B, Bu2) ≈ B*Bu2
        @test mul!(similar(B, T), Bu2, B) ≈ Bu2*B
    end
end

@testset "indexing partly initialized matrices" begin
    M = Matrix{BigFloat}(undef, 2, 2)
    U = UpperTriangular(M)
    @test iszero(U[2,1])
    L = LowerTriangular(M)
    @test iszero(L[1,2])
end

@testset "special printing of Lower/UpperTriangular" begin
    @test occursin(r"3×3 (LinearAlgebra\.)?LowerTriangular{Int64, Matrix{Int64}}:\n 2  ⋅  ⋅\n 2  2  ⋅\n 2  2  2",
                   sprint(show, MIME"text/plain"(), LowerTriangular(2ones(Int64,3,3))))
    @test occursin(r"3×3 (LinearAlgebra\.)?UnitLowerTriangular{Int64, Matrix{Int64}}:\n 1  ⋅  ⋅\n 2  1  ⋅\n 2  2  1",
                   sprint(show, MIME"text/plain"(), UnitLowerTriangular(2ones(Int64,3,3))))
    @test occursin(r"3×3 (LinearAlgebra\.)?UpperTriangular{Int64, Matrix{Int64}}:\n 2  2  2\n ⋅  2  2\n ⋅  ⋅  2",
                   sprint(show, MIME"text/plain"(), UpperTriangular(2ones(Int64,3,3))))
    @test occursin(r"3×3 (LinearAlgebra\.)?UnitUpperTriangular{Int64, Matrix{Int64}}:\n 1  2  2\n ⋅  1  2\n ⋅  ⋅  1",
                   sprint(show, MIME"text/plain"(), UnitUpperTriangular(2ones(Int64,3,3))))

    # don't access non-structural elements while displaying
    M = Matrix{BigFloat}(undef, 2, 2)
    @test sprint(show, UpperTriangular(M)) == "BigFloat[#undef #undef; 0.0 #undef]"
    @test sprint(show, LowerTriangular(M)) == "BigFloat[#undef 0.0; #undef #undef]"
end

@testset "adjoint/transpose triangular/vector multiplication" begin
    for elty in (Float64, ComplexF64), trity in (UpperTriangular, LowerTriangular)
        A1 = trity(rand(elty, 1, 1))
        b1 = rand(elty, 1)
        A4 = trity(rand(elty, 4, 4))
        b4 = rand(elty, 4)
        @test A1 * b1' ≈ Matrix(A1) * b1'
        @test_throws DimensionMismatch A4 * b4'
        @test A1 * transpose(b1) ≈ Matrix(A1) * transpose(b1)
        @test_throws DimensionMismatch A4 * transpose(b4)
        @test A1' * b1' ≈ Matrix(A1') * b1'
        @test_throws DimensionMismatch A4' * b4'
        @test A1' * transpose(b1) ≈  Matrix(A1') * transpose(b1)
        @test_throws DimensionMismatch A4' * transpose(b4)
        @test transpose(A1) * transpose(b1) ≈  Matrix(transpose(A1)) * transpose(b1)
        @test_throws DimensionMismatch transpose(A4) * transpose(b4)
        @test transpose(A1) * b1' ≈ Matrix(transpose(A1)) * b1'
        @test_throws DimensionMismatch transpose(A4) * b4'
        @test b1' * transpose(A1) ≈ b1' * Matrix(transpose(A1))
        @test b4' * transpose(A4) ≈ b4' * Matrix(transpose(A4))
        @test transpose(b1) * A1' ≈ transpose(b1) * Matrix(A1')
        @test transpose(b4) * A4' ≈ transpose(b4) * Matrix(A4')
    end
end

@testset "Error condition for powm" begin
    A = UpperTriangular(rand(ComplexF64, 10, 10))
    @test_throws ArgumentError LinearAlgebra.powm!(A, 2.2)
    A = LowerTriangular(rand(ComplexF64, 10, 10))
    At = copy(transpose(A))
    p = rand()
    @test LinearAlgebra.powm(A, p) == transpose(LinearAlgebra.powm!(At, p))
    @test_throws ArgumentError LinearAlgebra.powm(A, 2.2)
end

# Issue 35058
let A = [0.9999999999999998 4.649058915617843e-16 -1.3149405273715513e-16 9.9959579317056e-17; -8.326672684688674e-16 1.0000000000000004 2.9280733590254494e-16 -2.9993900031619594e-16; 9.43689570931383e-16 -1.339206523454095e-15 1.0000000000000007 -8.550505126287743e-16; -6.245004513516506e-16 -2.0122792321330962e-16 1.183061278035052e-16 1.0000000000000002],
    B = [0.09648289218436859 0.023497875751503007 0.0 0.0; 0.023497875751503007 0.045787575150300804 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    @test sqrt(A*B*A')^2 ≈ A*B*A'
end

@testset "one and oneunit for triangular" begin
    m = rand(4,4)
    function test_one_oneunit_triangular(a)
        b = Matrix(a)
        @test (@inferred a^1) == b^1
        @test (@inferred a^-1) ≈ b^-1
        @test one(a) == one(b)
        @test one(a)*a == a
        @test a*one(a) == a
        @test oneunit(a) == oneunit(b)
        @test oneunit(a) isa typeof(a)
    end
    for T in [UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular]
        a = T(m)
        test_one_oneunit_triangular(a)
    end
    # more complicated examples
    b = UpperTriangular(LowerTriangular(m))
    test_one_oneunit_triangular(b)
    c = UpperTriangular(Diagonal(rand(2)))
    test_one_oneunit_triangular(c)
end

@testset "LowerTriangular(Diagonal(...)) and friends (issue #28869)" begin
    for elty in (Float32, Float64, BigFloat, ComplexF32, ComplexF64, Complex{BigFloat}, Int)
        V = elty ≡ Int ? rand(1:10, 5) : elty.(randn(5))
        D = Diagonal(V)
        for dty in (UpperTriangular, LowerTriangular)
            A = dty(D)
            @test A * A' == D * D'
        end
    end
end

@testset "tril!/triu! for non-bitstype matrices" begin
    @testset "numeric" begin
        M = Matrix{BigFloat}(undef, 3, 3)
        tril!(M)
        L = LowerTriangular(ones(3,3))
        copytrito!(M, L, 'L')
        @test M == L

        M = Matrix{BigFloat}(undef, 3, 3)
        triu!(M)
        U = UpperTriangular(ones(3,3))
        copytrito!(M, U, 'U')
        @test M == U
    end
    @testset "array elements" begin
        M = fill(ones(2,2), 4, 4)
        tril!(M)
        L = LowerTriangular(fill(fill(2,2,2),4,4))
        copytrito!(M, L, 'L')
        @test M == L

        M = fill(ones(2,2), 4, 4)
        triu!(M)
        U = UpperTriangular(fill(fill(2,2,2),4,4))
        copytrito!(M, U, 'U')
        @test M == U
    end
end

@testset "avoid matmul ambiguities with ::MyMatrix * ::AbstractMatrix" begin
    A = [i+j for i in 1:2, j in 1:2]
    S = SizedArrays.SizedArray{(2,2)}(A)
    U = UpperTriangular(ones(2,2))
    @test S * U == A * U
    @test U * S == U * A
    C1, C2 = zeros(2,2), zeros(2,2)
    @test mul!(C1, S, U) == mul!(C2, A, U)
    @test mul!(C1, S, U, 1, 2) == mul!(C2, A, U, 1 ,2)
    @test mul!(C1, U, S) == mul!(C2, U, A)
    @test mul!(C1, U, S, 1, 2) == mul!(C2, U, A, 1 ,2)

    v = [i for i in 1:2]
    sv = SizedArrays.SizedArray{(2,)}(v)
    @test U * sv == U * v
    C1, C2 = zeros(2), zeros(2)
    @test mul!(C1, U, sv) == mul!(C2, U, v)
    @test mul!(C1, U, sv, 1, 2) == mul!(C2, U, v, 1 ,2)
end

@testset "custom axes" begin
    SZA = SizedArrays.SizedArray{(2,2)}([1 2; 3 4])
    for T in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
        S = T(SZA)
        r = SizedArrays.SOneTo(2)
        @test axes(S) === (r,r)
    end
end

@testset "immutable and non-strided parent" begin
    F = FillArrays.Fill(2, (4,4))
    for UT in (UnitUpperTriangular, UnitLowerTriangular)
        U = UT(F)
        @test -U == -Array(U)
    end

    F = FillArrays.Fill(3im, (4,4))
    for U in (UnitUpperTriangular(F), UnitLowerTriangular(F))
        @test imag(F) == imag(collect(F))
    end

    @testset "copyto!" begin
        for T in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
            @test Matrix(T(F)) == T(F)
        end
        @test copyto!(zeros(eltype(F), length(F)), UpperTriangular(F)) == vec(UpperTriangular(F))
    end
end

@testset "error paths" begin
    A = zeros(1,1); B = zeros(2,2)
    @testset "inplace mul scaling with incompatible sizes" begin
        for T in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
            @test_throws DimensionMismatch mul!(T(A), T(B), 3)
            @test_throws DimensionMismatch mul!(T(A), 3, T(B))
        end
    end
    @testset "copyto with incompatible sizes" begin
        for T in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
            @test_throws BoundsError copyto!(T(A), T(B))
        end
    end
end

@testset "uppertriangular/lowertriangular" begin
    M = rand(2,2)
    @test LinearAlgebra.uppertriangular(M) === UpperTriangular(M)
    @test LinearAlgebra.lowertriangular(M) === LowerTriangular(M)
    @test LinearAlgebra.uppertriangular(UnitUpperTriangular(M)) === UnitUpperTriangular(M)
    @test LinearAlgebra.lowertriangular(UnitLowerTriangular(M)) === UnitLowerTriangular(M)
end

@testset "arithmetic with partly uninitialized matrices" begin
    @testset "$(typeof(A))" for A in (Matrix{BigFloat}(undef,2,2), Matrix{Complex{BigFloat}}(undef,2,2)')
        A[2,1] = eltype(A) <: Complex ? 4 + 3im : 4
        B = Matrix{eltype(A)}(undef, size(A))
        for MT in (LowerTriangular, UnitLowerTriangular)
            if MT == LowerTriangular
                A[1,1] = A[2,2] = eltype(A) <: Complex ? 4 + 3im : 4
            end
            L = MT(A)
            B .= 0
            copyto!(B, L)
            @test copy(L) == B
            @test L * 2 == 2 * L == 2B
            @test L/2 == B/2
            @test 2\L == 2\B
            @test real(L) == real(B)
            @test imag(L) == imag(B)
            if MT == LowerTriangular
                @test isa(kron(L,L), MT)
            end
            @test kron(L,L) == kron(B,B)
            @test transpose!(MT(copy(A))) == transpose(L) broken=!(A isa Matrix)
            @test adjoint!(MT(copy(A))) == adjoint(L) broken=!(A isa Matrix)
        end
    end

    @testset "$(typeof(A))" for A in (Matrix{BigFloat}(undef,2,2), Matrix{Complex{BigFloat}}(undef,2,2)')
        A[1,2] = eltype(A) <: Complex ? 4 + 3im : 4
        B = Matrix{eltype(A)}(undef, size(A))
        for MT in (UpperTriangular, UnitUpperTriangular)
            if MT == UpperTriangular
                A[1,1] = A[2,2] = eltype(A) <: Complex ? 4 + 3im : 4
            end
            U = MT(A)
            B .= 0
            copyto!(B, U)
            @test copy(U) == B
            @test U * 2 == 2 * U == 2B
            @test U/2 == B/2
            @test 2\U == 2\B
            @test real(U) == real(B)
            @test imag(U) == imag(B)
            if MT == UpperTriangular
                @test isa(kron(U,U), MT)
            end
            @test kron(U,U) == kron(B,B)
            @test transpose!(MT(copy(A))) == transpose(U) broken=!(A isa Matrix)
            @test adjoint!(MT(copy(A))) == adjoint(U) broken=!(A isa Matrix)
        end
    end
end

@testset "kron with triangular matrices of matrices" begin
    for T in (UpperTriangular, LowerTriangular)
        t = T(fill(ones(2,2), 2, 2))
        m = Matrix(t)
        @test isa(kron(t,t), T)
        @test kron(t, t) ≈ kron(m, m)
    end
end

@testset "kron with triangular matrices of mixed eltypes" begin
    for T in (UpperTriangular, LowerTriangular)
        U = T(Matrix{Union{Missing,Int}}(fill(2, 2, 2)))
        U[1, 1] = missing
        @test kron(U, U)[2, 3] == 0
        @test kron(U, U)[3, 2] == 0
    end
end

@testset "copyto! tests" begin
    @testset "copyto! with aliasing (#39460)" begin
        M = Matrix(reshape(1:36, 6, 6))
        @testset for T in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
            A = T(view(M, 1:5, 1:5))
            A2 = copy(A)
            B = T(view(M, 2:6, 2:6))
            @test copyto!(B, A) == A2
        end
    end

    @testset "copyto! with different matrix types" begin
        M1 = Matrix(reshape(1:36, 6, 6))
        M2 = similar(M1)
        # these copies always work
        @testset for (Tdest, Tsrc) in (
                            (UpperTriangular, UnitUpperTriangular),
                            (UpperTriangular, UpperTriangular),
                            (LowerTriangular, UnitLowerTriangular),
                            (LowerTriangular, LowerTriangular),
                            (UnitUpperTriangular, UnitUpperTriangular),
                            (UnitLowerTriangular, UnitLowerTriangular)
                        )

            M2 .= 0
            copyto!(Tdest(M2), Tsrc(M1))
            @test Tdest(M2) == Tsrc(M1)
        end
        # these copies only work if the source has a unit diagonal
        M3 = copy(M1)
        M3[diagind(M3)] .= 1
        @testset for (Tdest, Tsrc) in (
                            (UnitUpperTriangular, UpperTriangular),
                            (UnitLowerTriangular, LowerTriangular),
                        )

            M2 .= 0
            copyto!(Tdest(M2), Tsrc(M3))
            @test Tdest(M2) == Tsrc(M3)
            @test_throws ArgumentError copyto!(Tdest(M2), Tsrc(M1))
        end
        # these copies work even when the parent of the source isn't initialized along the diagonal
        @testset for (T, TU) in ((UpperTriangular, UnitUpperTriangular),
                                    (LowerTriangular, UnitLowerTriangular))
            M1 = Matrix{BigFloat}(undef, 3, 3)
            M2 = similar(M1)
            if TU == UnitUpperTriangular
                M1[1,2] = M1[1,3] = M1[2,3] = 2
            else
                M1[2,1] = M1[3,1] = M1[3,2] = 2
            end
            for TD in (T, TU)
                M2 .= 0
                copyto!(T(M2), TU(M1))
                @test T(M2) == TU(M1)
            end
        end
    end

    @testset "copyto! with different sizes" begin
        Ap = zeros(3,3)
        Bp = rand(2,2)
        @testset for T in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
            A = T(Ap)
            B = T(Bp)
            @test_throws ArgumentError copyto!(A, B)
        end
        @testset "error message" begin
            A = UpperTriangular(Ap)
            B = UpperTriangular(Bp)
            @test_throws "cannot set index in the lower triangular part" copyto!(A, B)

            A = LowerTriangular(Ap)
            B = LowerTriangular(Bp)
            @test_throws "cannot set index in the upper triangular part" copyto!(A, B)
        end
    end
end

@testset "getindex with Integers" begin
    M = reshape(1:4,2,2)
    for Ttype in (UpperTriangular, UnitUpperTriangular)
        T = Ttype(M)
        @test_throws "invalid index" T[2, true]
        @test T[1,2] == T[Int8(1),UInt16(2)] == T[big(1), Int16(2)]
    end
    for Ttype in (LowerTriangular, UnitLowerTriangular)
        T = Ttype(M)
        @test_throws "invalid index" T[true, 2]
        @test T[2,1] == T[Int8(2),UInt16(1)] == T[big(2), Int16(1)]
    end
end

@testset "type-stable eigvecs" begin
    D = Float64[1 0; 0 2]
    V = @inferred eigvecs(UpperTriangular(D))
    @test V == Diagonal([1, 1])
end

@testset "preserve structure in scaling by NaN" begin
    M = rand(Int8,2,2)
    for (Ts, TD) in (((UpperTriangular, UnitUpperTriangular), UpperTriangular),
                    ((LowerTriangular, UnitLowerTriangular), LowerTriangular))
        for T in Ts
            U = T(M)
            for V in (U * NaN, NaN * U, U / NaN, NaN \ U)
                @test V isa TD{Float64, Matrix{Float64}}
                @test all(isnan, diag(V))
            end
        end
    end
end

@testset "eigvecs for AbstractTriangular" begin
    S = SizedArrays.SizedArray{(3,3)}(reshape(1:9,3,3))
    for T in (UpperTriangular, UnitUpperTriangular,
                LowerTriangular, UnitLowerTriangular)
        U = T(S)
        V = eigvecs(U)
        λ = eigvals(U)
        @test U * V ≈ V * Diagonal(λ)

        MU = MyTriangular(U)
        V = eigvecs(U)
        λ = eigvals(U)
        @test MU * V ≈ V * Diagonal(λ)
    end
end

@testset "(l/r)mul! and (l/r)div! for generic triangular" begin
    @testset for T in (UpperTriangular, LowerTriangular, UnitUpperTriangular, UnitLowerTriangular)
        M = MyTriangular(T(rand(4,4)))
        D = Diagonal(randn(4))
        A = rand(4,4)
        Ac = similar(A)
        @testset "lmul!" begin
            Ac .= A
            lmul!(M, Ac)
            @test Ac ≈ M * A
        end
        @testset "rmul!" begin
            Ac .= A
            rmul!(Ac, M)
            @test Ac ≈ A * M
        end
        @testset "ldiv!" begin
            Ac .= A
            ldiv!(M, Ac)
            @test Ac ≈ M \ A
        end
        @testset "rdiv!" begin
            Ac .= A
            rdiv!(Ac, M)
            @test Ac ≈ A / M
        end
        @testset "diagonal mul" begin
            @test D * M ≈ D * M.data
            @test M * D ≈ M.data * D
        end
    end
end

@testset "istriu/istril forwards to parent" begin
    @testset "$(nameof(typeof(M)))" for M in [Tridiagonal(rand(n-1), rand(n), rand(n-1)),
                Tridiagonal(zeros(n-1), zeros(n), zeros(n-1)),
                Diagonal(randn(n)),
                Diagonal(zeros(n)),
                ]
        @testset for TriT in (UpperTriangular, UnitUpperTriangular, LowerTriangular, UnitLowerTriangular)
            U = TriT(M)
            A = Array(U)
            for k in -n:n
                @test istriu(U, k) == istriu(A, k)
                @test istril(U, k) == istril(A, k)
            end
        end
    end
    z = zeros(n,n)
    @testset for TriT in (UpperTriangular, UnitUpperTriangular, LowerTriangular, UnitLowerTriangular)
        P = Matrix{BigFloat}(undef, n, n)
        copytrito!(P, z, TriT <: Union{UpperTriangular, UnitUpperTriangular} ? 'U' : 'L')
        U = TriT(P)
        A = Array(U)
        @testset for k in -n:n
            @test istriu(U, k) == istriu(A, k)
            @test istril(U, k) == istril(A, k)
        end
    end

    @testset "Union eltype" begin
        M = Matrix{Union{Int,Missing}}(missing,2,2)
        U = triu(M)
        @test iszero(U[2,1])
        U = tril(M)
        @test iszero(U[1,2])
    end
end

@testset "indexing with a BandIndex" begin
    # these tests should succeed even if the linear index along
    # the band isn't a constant, or type-inferred at all
    M = rand(Int,2,2)
    f(A,j, v::Val{n}) where {n} = Val(A[BandIndex(n,j)])
    function common_tests(M, ind)
        j = ind[]
        @test @inferred(f(UpperTriangular(M), j, Val(-1))) == Val(0)
        @test @inferred(f(UnitUpperTriangular(M), j, Val(-1))) == Val(0)
        @test @inferred(f(UnitUpperTriangular(M), j, Val(0))) == Val(1)
        @test @inferred(f(LowerTriangular(M), j, Val(1))) == Val(0)
        @test @inferred(f(UnitLowerTriangular(M), j, Val(1))) == Val(0)
        @test @inferred(f(UnitLowerTriangular(M), j, Val(0))) == Val(1)
    end
    common_tests(M, Any[1])

    M = Diagonal([1,2])
    common_tests(M, Any[1])
    # extra tests for banded structure of the parent
    for T in (UpperTriangular, UnitUpperTriangular)
        @test @inferred(f(T(M), 1, Val(1))) == Val(0)
    end
    for T in (LowerTriangular, UnitLowerTriangular)
        @test @inferred(f(T(M), 1, Val(-1))) == Val(0)
    end

    M = Tridiagonal([1,2], [1,2,3], [1,2])
    common_tests(M, Any[1])
    for T in (UpperTriangular, UnitUpperTriangular)
        @test @inferred(f(T(M), 1, Val(2))) == Val(0)
    end
    for T in (LowerTriangular, UnitLowerTriangular)
        @test @inferred(f(T(M), 1, Val(-2))) == Val(0)
    end
end

@testset "indexing uses diagzero" begin
    @testset "block matrix" begin
        M = reshape([zeros(2,2), zeros(4,2), zeros(2,3), zeros(4,3)],2,2)
        U = UpperTriangular(M)
        @test [size(x) for x in U] == [size(x) for x in M]
    end
    @testset "Union eltype" begin
        M = Matrix{Union{Int,Missing}}(missing,4,4)
        U = UpperTriangular(M)
        @test iszero(U[3,1])
    end
end

@testset "addition/subtraction of mixed triangular" begin
    for A in (Hermitian(rand(4, 4)), Diagonal(rand(5)))
        for T in (UpperTriangular, LowerTriangular,
                UnitUpperTriangular, UnitLowerTriangular)
            B = T(A)
            M = Matrix(B)
            R = B - B'
            if A isa Diagonal
                @test R isa Diagonal
            end
            @test R == M - M'
            R = B + B'
            if A isa Diagonal
                @test R isa Diagonal
            end
            @test R == M + M'
            C = MyTriangular(B)
            @test C - C' == M - M'
            @test C + C' == M + M'
        end
    end
    @testset "unfilled parent" begin
        @testset for T in (UpperTriangular, LowerTriangular,
                UnitUpperTriangular, UnitLowerTriangular)
            F = Matrix{BigFloat}(undef, 2, 2)
            B = T(F)
            isupper = B isa Union{UpperTriangular, UnitUpperTriangular}
            B[1+!isupper, 1+isupper] = 2
            if !(B isa Union{UnitUpperTriangular, UnitLowerTriangular})
                B[1,1] = B[2,2] = 3
            end
            M = Matrix(B)
            @test B - B' == M - M'
            @test B + B' == M + M'
            @test B - copy(B') == M - M'
            @test B + copy(B') == M + M'
            C = MyTriangular(B)
            @test C - C' == M - M'
            @test C + C' == M + M'
        end
    end
end

@testset "log_quasitriu with internal scaling s=0 (issue #54833)" begin
    M = [0.9949357359852791 -0.015567763143324862 -0.09091193493947397 -0.03994428739762443 0.07338356301650806;
    0.011813655598647289 0.9968988574699793 -0.06204555000202496 0.04694097614450692 0.09028834462782365;
    0.092737943594701 0.059546719185135925 0.9935850721633324 0.025348893985651405 -0.018530261590167685;
    0.0369187299165628 -0.04903571106913449 -0.025962938675946543 0.9977767446862031 0.12901494726320517;
    0.0 0.0 0.0 0.0 1.0]

    @test exp(log(M)) ≈ M
end

@testset "copytrito!" begin
    for T in (UpperTriangular, LowerTriangular)
        M = Matrix{BigFloat}(undef, 2, 2)
        M[1,1] = M[2,2] = 3
        U = T(M)
        isupper = U isa UpperTriangular
        M[1+!isupper, 1+isupper] = 4
        uplo, loup = U isa UpperTriangular ? ('U', 'L') : ('L', 'U' )
        @test copytrito!(similar(U), U, uplo) == U
        @test copytrito!(zero(M), U, uplo) == U
        @test copytrito!(similar(U), Array(U), uplo) == U
        @test copytrito!(zero(U), U, loup) == Diagonal(U)
        @test copytrito!(similar(U), MyTriangular(U), uplo) == U
        @test copytrito!(zero(M), MyTriangular(U), uplo) == U
        Ubig = T(similar(M, (3,3)))
        copytrito!(Ubig, U, uplo)
        @test Ubig[axes(U)...] == U
    end
end

@testset "(l/r)mul! and (l/r)div! for non-contiguous arrays" begin
    U = UpperTriangular(reshape(collect(3:27.0),5,5))
    b = float.(1:10)
    b2 = copy(b); b2v = view(b2, 1:2:9); b2vc = copy(b2v)
    @test lmul!(U, b2v) == lmul!(U, b2vc)
    b2 = copy(b); b2v = view(b2, 1:2:9); b2vc = copy(b2v)
    @test ldiv!(U, b2v) ≈ ldiv!(U, b2vc)
    B = float.(collect(reshape(1:100, 10,10)))
    B2 = copy(B); B2v = view(B2, 1:2:9, 1:5); B2vc = copy(B2v)
    @test lmul!(U, B2v) == lmul!(U, B2vc)
    B2 = copy(B); B2v = view(B2, 1:2:9, 1:5); B2vc = copy(B2v)
    @test rmul!(B2v, U) == rmul!(B2vc, U)
    B2 = copy(B); B2v = view(B2, 1:2:9, 1:5); B2vc = copy(B2v)
    @test ldiv!(U, B2v) ≈ ldiv!(U, B2vc)
    B2 = copy(B); B2v = view(B2, 1:2:9, 1:5); B2vc = copy(B2v)
    @test rdiv!(B2v, U) ≈ rdiv!(B2vc, U)
end

@testset "error messages in matmul with mismatched matrix sizes" begin
    for T in (Int, Float64)
        A = UpperTriangular(ones(T,2,2))
        B = ones(T,3,3)
        C = similar(B)
        @test_throws "incompatible dimensions for matrix multiplication" mul!(C, A, B)
        @test_throws "incompatible dimensions for matrix multiplication" mul!(C, B, A)
        B = Array(A)
        C = similar(B, (4,4))
        @test_throws "incompatible destination size" mul!(C, A, B)
        @test_throws "incompatible destination size" mul!(C, B, A)
    end
end

end # module TestTriangular

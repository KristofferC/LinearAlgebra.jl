# This file is a part of Julia. License is MIT: https://julialang.org/license

# The following test block tries to call all methods in base/linalg/triangular.jl in order for a combination of input element types. Keep the ordering when adding code.
function test_triangular(elty1_types)
    n = 9
    @testset for elty1 in elty1_types
        # Begin loop for first Triangular matrix
        @testset for (t1, uplo1) in ((UpperTriangular, :U),
            (UnitUpperTriangular, :U),
            (LowerTriangular, :L),
            (UnitLowerTriangular, :L))

            # Construct test matrix
            A1 = t1(elty1 == Int ? rand(1:7, n, n) : convert(Matrix{elty1}, (elty1 <: Complex ? complex.(randn(n, n), randn(n, n)) : randn(n, n)) |> t -> cholesky(t't).U |> t -> uplo1 === :U ? t : copy(t')))
            M1 = Matrix(A1)
            @test t1(A1) === A1
            @test t1{elty1}(A1) === A1
            # test the ctor works for AbstractMatrix
            symm = Symmetric(rand(Int8, n, n))
            t1s = t1{elty1}(symm)
            @test typeof(t1s) == t1{elty1,Symmetric{elty1,Matrix{elty1}}}
            t1t = t1{elty1}(t1(rand(Int8, n, n)))
            @test typeof(t1t) == t1{elty1,Matrix{elty1}}

            # Convert
            @test convert(AbstractMatrix{elty1}, A1) == A1
            @test convert(Matrix, A1) == A1
            @test t1{elty1}(convert(AbstractMatrix{elty1}, A1)) == A1

            # full!
            @test full!(copy(A1)) == A1

            # similar
            @test isa(similar(A1), t1)
            @test eltype(similar(A1)) == elty1
            @test isa(similar(A1, Int), t1)
            @test eltype(similar(A1, Int)) == Int
            @test isa(similar(A1, (3, 2)), Matrix{elty1})
            @test isa(similar(A1, Int, (3, 2)), Matrix{Int})

            #copyto!
            simA1 = similar(A1)
            copyto!(simA1, A1)
            @test simA1 == A1

            # getindex
            let mA1 = M1
                # linear indexing
                for i in 1:length(A1)
                    @test A1[i] == mA1[i]
                end
                # cartesian indexing
                for i in 1:size(A1, 1), j in 1:size(A1, 2)
                    @test A1[i, j] == mA1[i, j]
                end
            end
            @test isa(A1[2:4, 1], Vector)


            # setindex! (and copy)
            A1c = copy(A1)
            for i = 1:size(A1, 1)
                for j = 1:size(A1, 2)
                    if uplo1 === :U
                        if i > j
                            A1c[i, j] = 0
                            @test_throws ArgumentError A1c[i, j] = 1
                        elseif i == j && t1 == UnitUpperTriangular
                            A1c[i, j] = 1
                            @test_throws ArgumentError A1c[i, j] = 0
                        else
                            A1c[i, j] = 0
                            @test A1c[i, j] == 0
                        end
                    else
                        if i < j
                            A1c[i, j] = 0
                            @test_throws ArgumentError A1c[i, j] = 1
                        elseif i == j && t1 == UnitLowerTriangular
                            A1c[i, j] = 1
                            @test_throws ArgumentError A1c[i, j] = 0
                        else
                            A1c[i, j] = 0
                            @test A1c[i, j] == 0
                        end
                    end
                end
            end

            # istril/istriu
            if uplo1 === :L
                @test istril(A1)
                @test !istriu(A1)
                @test istriu(A1')
                @test istriu(transpose(A1))
                @test !istril(A1')
                @test !istril(transpose(A1))
            else
                @test istriu(A1)
                @test !istril(A1)
                @test istril(A1')
                @test istril(transpose(A1))
                @test !istriu(A1')
                @test !istriu(transpose(A1))
            end
            M = copy(parent(A1))
            for trans in (adjoint, transpose), k in -1:1
                triu!(M, k)
                @test istril(trans(M), -k) == istril(copy(trans(M)), -k) == true
            end
            M = copy(parent(A1))
            for trans in (adjoint, transpose), k in 1:-1:-1
                tril!(M, k)
                @test istriu(trans(M), -k) == istriu(copy(trans(M)), -k) == true
            end

            #tril/triu
            if uplo1 === :L
                @test tril(A1, 0) == A1
                @test tril(A1, -1) == LowerTriangular(tril(M1, -1))
                @test tril(A1, 1) == t1(tril(tril(M1, 1)))
                @test tril(A1, -n - 2) == zeros(size(A1))
                @test tril(A1, n) == A1
                @test triu(A1, 0) == t1(diagm(0 => diag(A1)))
                @test triu(A1, -1) == t1(tril(triu(A1.data, -1)))
                @test triu(A1, 1) == zeros(size(A1)) # or just @test iszero(triu(A1,1))?
                @test triu(A1, -n) == A1
                @test triu(A1, n + 2) == zeros(size(A1))
            else
                @test triu(A1, 0) == A1
                @test triu(A1, 1) == UpperTriangular(triu(M1, 1))
                @test triu(A1, -1) == t1(triu(triu(M1, -1)))
                @test triu(A1, -n) == A1
                @test triu(A1, n + 2) == zeros(size(A1))
                @test tril(A1, 0) == t1(diagm(0 => diag(A1)))
                @test tril(A1, 1) == t1(triu(tril(A1.data, 1)))
                @test tril(A1, -1) == zeros(size(A1)) # or just @test iszero(tril(A1,-1))?
                @test tril(A1, -n - 2) == zeros(size(A1))
                @test tril(A1, n) == A1
            end

            # factorize
            @test factorize(A1) == A1

            # [c]transpose[!] (test views as well, see issue #14317)
            let vrange = 1:n-1, viewA1 = t1(view(A1.data, vrange, vrange))
                # transpose
                @test copy(transpose(A1)) == transpose(M1)
                @test copy(transpose(viewA1)) == transpose(Matrix(viewA1))
                # adjoint
                @test copy(A1') == M1'
                @test copy(viewA1') == Matrix(viewA1)'
                # transpose!
                @test transpose!(copy(A1)) == transpose(A1)
                @test typeof(transpose!(copy(A1))).name == typeof(transpose(A1)).name
                @test transpose!(t1(view(copy(A1).data, vrange, vrange))) == transpose(viewA1)
                # adjoint!
                @test adjoint!(copy(A1)) == adjoint(A1)
                @test typeof(adjoint!(copy(A1))).name == typeof(adjoint(A1)).name
                @test adjoint!(t1(view(copy(A1).data, vrange, vrange))) == adjoint(viewA1)
            end

            # diag
            @test diag(A1) == diag(M1)

            # tr
            @test tr(A1)::elty1 == tr(M1)

            # real
            @test real(A1) == real(M1)
            @test imag(A1) == imag(M1)
            @test abs.(A1) == abs.(M1)

            # zero
            if A1 isa UpperTriangular || A1 isa LowerTriangular
                @test zero(A1) == zero(parent(A1))
            end

            # Unary operations
            @test -A1 == -M1

            # copy and copyto! (test views as well, see issue #14317)
            let vrange = 1:n-1, viewA1 = t1(view(A1.data, vrange, vrange))
                # copy
                @test copy(A1) == copy(M1)
                @test copy(viewA1) == copy(Matrix(viewA1))
                # copyto!
                B = similar(A1)
                copyto!(B, A1)
                @test B == A1
                B = similar(copy(transpose(A1)))
                copyto!(B, copy(transpose(A1)))
                @test B == copy(transpose(A1))
                B = similar(viewA1)
                copyto!(B, viewA1)
                @test B == viewA1
                B = similar(copy(transpose(viewA1)))
                copyto!(B, copy(transpose(viewA1)))
                @test B == transpose(viewA1)
            end

            #exp/log
            if elty1 ∈ (Float32, Float64, ComplexF32, ComplexF64)
                @test exp(Matrix(log(A1))) ≈ A1
            end

            # scale
            if (t1 == UpperTriangular || t1 == LowerTriangular)
                unitt = istriu(A1) ? UnitUpperTriangular : UnitLowerTriangular
                if elty1 == Int
                    cr = 2
                else
                    cr = 0.5
                end
                ci = cr * im
                if elty1 <: Real
                    A1tmp = copy(A1)
                    rmul!(A1tmp, cr)
                    @test A1tmp == cr * A1
                    A1tmp = copy(A1)
                    lmul!(cr, A1tmp)
                    @test A1tmp == cr * A1
                    A1tmp = copy(A1)
                    A2tmp = unitt(A1)
                    mul!(A1tmp, A2tmp, cr)
                    @test A1tmp == cr * A2tmp
                    A1tmp = copy(A1)
                    A2tmp = unitt(A1)
                    mul!(A1tmp, cr, A2tmp)
                    @test A1tmp == cr * A2tmp

                    A1tmp .= A1
                    @test mul!(A1tmp, A2tmp, cr, 0, 2) == 2A1
                    A1tmp .= A1
                    @test mul!(A1tmp, cr, A2tmp, 0, 2) == 2A1
                else
                    A1tmp = copy(A1)
                    rmul!(A1tmp, ci)
                    @test A1tmp == ci * A1
                    A1tmp = copy(A1)
                    lmul!(ci, A1tmp)
                    @test A1tmp == ci * A1
                    A1tmp = copy(A1)
                    A2tmp = unitt(A1)
                    mul!(A1tmp, ci, A2tmp)
                    @test A1tmp == ci * A2tmp
                    A1tmp = copy(A1)
                    A2tmp = unitt(A1)
                    mul!(A1tmp, A2tmp, ci)
                    @test A1tmp == A2tmp * ci
                end
            end

            # generalized dot
            for eltyb in (Float32, Float64, BigFloat, ComplexF32, ComplexF64, Complex{BigFloat})
                b1 = convert(Vector{eltyb}, (elty1 <: Complex ? real(A1) : A1) * fill(1.0, n))
                b2 = convert(Vector{eltyb}, (elty1 <: Complex ? real(A1) : A1) * randn(n))
                @test dot(b1, A1, b2) ≈ dot(A1'b1, b2) atol = sqrt(max(eps(real(float(one(elty1)))), eps(real(float(one(eltyb)))))) * n * n
            end

            # Binary operations
            @test A1 * 0.5 == M1 * 0.5
            @test 0.5 * A1 == 0.5 * M1
            @test A1 / 0.5 == M1 / 0.5
            @test 0.5 \ A1 == 0.5 \ M1

            # inversion
            @test inv(A1) ≈ inv(lu(M1))
            inv(M1) # issue #11298
            @test isa(inv(A1), t1)
            # make sure the call to LAPACK works right
            if elty1 <: BlasFloat
                @test LinearAlgebra.inv!(copy(A1)) ≈ inv(lu(M1))
            end

            # Determinant
            @test det(A1) ≈ det(lu(M1)) atol = sqrt(eps(real(float(one(elty1))))) * n * n
            @test logdet(A1) ≈ logdet(lu(M1)) atol = sqrt(eps(real(float(one(elty1))))) * n * n
            lada, ladb = logabsdet(A1)
            flada, fladb = logabsdet(lu(M1))
            @test lada ≈ flada atol = sqrt(eps(real(float(one(elty1))))) * n * n
            @test ladb ≈ fladb atol = sqrt(eps(real(float(one(elty1))))) * n * n

            # Matrix square root
            @test sqrt(A1) |> (t -> (t * t)::typeof(t)) ≈ A1

            # naivesub errors
            @test_throws DimensionMismatch ldiv!(A1, Vector{elty1}(undef, n + 1))

            # eigenproblems
            if !(elty1 in (BigFloat, Complex{BigFloat})) # Not handled yet
                vals, vecs = eigen(A1)
                if (t1 == UpperTriangular || t1 == LowerTriangular) && elty1 != Int # Cannot really handle degenerate eigen space and Int matrices will probably have repeated eigenvalues.
                    @test vecs * diagm(0 => vals) / vecs ≈ A1 atol = sqrt(eps(float(real(one(vals[1]))))) * (opnorm(A1, Inf) * n)^2
                end
            end

            # Condition number tests - can be VERY approximate
            if elty1 <: BlasFloat
                for p in (1.0, Inf)
                    @test cond(A1, p) ≈ cond(A1, p) atol = (cond(A1, p) + cond(A1, p))
                end
                @test cond(A1, 2) == cond(M1, 2)
            end

            if !(elty1 in (BigFloat, Complex{BigFloat})) # Not implemented yet
                svd(A1)
                elty1 <: BlasFloat && svd!(copy(A1))
                svdvals(A1)
            end

            @test ((A1 * A1)::t1) ≈ M1 * M1
            @test ((A1 / A1)::t1) ≈ M1 / M1
            @test ((A1 \ A1)::t1) ≈ M1 \ M1

            # Begin loop for second Triangular matrix
            @testset for elty2 in (Float32, Float64, BigFloat, ComplexF32, ComplexF64, Complex{BigFloat}, Int)
                # Only test methods for the same element type and a single combination of mixed element types
                # to avoid too much compilation
                if !(elty1 == elty2 || elty1 ∈ (ComplexF32, Int) || elty1 ∈ (ComplexF32, Int))
                    continue
                end
                @testset for (t2, uplo2) in ((UpperTriangular, :U),
                    (UnitUpperTriangular, :U),
                    (LowerTriangular, :L),
                    (UnitLowerTriangular, :L))

                    A2 = t2(elty2 == Int ? rand(1:7, n, n) : convert(Matrix{elty2}, (elty2 <: Complex ? complex.(randn(n, n), randn(n, n)) : randn(n, n)) |> t -> cholesky(t't).U |> t -> uplo2 === :U ? t : copy(t')))
                    M2 = Matrix(A2)
                    # Convert
                    if elty1 <: Real && !(elty2 <: Integer)
                        @test convert(AbstractMatrix{elty2}, A1) == t1(convert(Matrix{elty2}, A1.data))
                    elseif elty2 <: Real && !(elty1 <: Integer)
                        @test_throws InexactError convert(AbstractMatrix{elty2}, A1) == t1(convert(Matrix{elty2}, A1.data))
                    end

                    # Binary operations
                    @test A1 + A2 == M1 + M2
                    @test A1 - A2 == M1 - M2
                    @test kron(A1, A2) == kron(M1, M2)

                    # Triangular-Triangular multiplication and division
                    @test A1 * A2 ≈ M1 * M2
                    @test transpose(A1) * A2 ≈ transpose(M1) * M2
                    @test transpose(A1) * adjoint(A2) ≈ transpose(M1) * adjoint(M2)
                    @test adjoint(A1) * transpose(A2) ≈ adjoint(M1) * transpose(M2)
                    @test A1'A2 ≈ M1'M2
                    @test A1 * transpose(A2) ≈ M1 * transpose(M2)
                    @test A1 * A2' ≈ M1 * M2'
                    @test transpose(A1) * transpose(A2) ≈ transpose(M1) * transpose(M2)
                    @test A1'A2' ≈ M1'M2'
                    @test A1 / A2 ≈ M1 / M2
                    @test A1 \ A2 ≈ M1 \ M2
                    if uplo1 === :U && uplo2 === :U
                        if t1 === UnitUpperTriangular && t2 === UnitUpperTriangular
                            @test A1 * A2 isa UnitUpperTriangular
                            @test A1 / A2 isa UnitUpperTriangular
                            elty1 == Int && elty2 == Int && @test eltype(A1 / A2) == Int
                            @test A1 \ A2 isa UnitUpperTriangular
                            elty1 == Int && elty2 == Int && @test eltype(A1 \ A2) == Int
                        else
                            @test A1 * A2 isa UpperTriangular
                            @test A1 / A2 isa UpperTriangular
                            elty1 == Int && elty2 == Int && t2 === UnitUpperTriangular && @test eltype(A1 / A2) == Int
                            @test A1 \ A2 isa UpperTriangular
                            elty1 == Int && elty2 == Int && t1 === UnitUpperTriangular && @test eltype(A1 \ A2) == Int
                        end
                    elseif uplo1 === :L && uplo2 === :L
                        if t1 === UnitLowerTriangular && t2 === UnitLowerTriangular
                            @test A1 * A2 isa UnitLowerTriangular
                            @test A1 / A2 isa UnitLowerTriangular
                            elty1 == Int && elty2 == Int && @test eltype(A1 / A2) == Int
                            @test A1 \ A2 isa UnitLowerTriangular
                            elty1 == Int && elty2 == Int && @test eltype(A1 \ A2) == Int
                        else
                            @test A1 * A2 isa LowerTriangular
                            @test A1 / A2 isa LowerTriangular
                            elty1 == Int && elty2 == Int && t2 === UnitLowerTriangular && @test eltype(A1 / A2) == Int
                            @test A1 \ A2 isa LowerTriangular
                            elty1 == Int && elty2 == Int && t1 === UnitLowerTriangular && @test eltype(A1 \ A2) == Int
                        end
                    end
                    offsizeA = Matrix{Float64}(I, n + 1, n + 1)
                    @test_throws DimensionMismatch offsizeA / A2
                    @test_throws DimensionMismatch offsizeA / transpose(A2)
                    @test_throws DimensionMismatch offsizeA / A2'
                    @test_throws DimensionMismatch offsizeA * A2
                    @test_throws DimensionMismatch offsizeA * transpose(A2)
                    @test_throws DimensionMismatch offsizeA * A2'
                    @test_throws DimensionMismatch transpose(A2) * offsizeA
                    @test_throws DimensionMismatch A2' * offsizeA
                    @test_throws DimensionMismatch A2 * offsizeA
                    if (uplo1 == uplo2 && elty1 == elty2 != Int && t1 != UnitLowerTriangular && t1 != UnitUpperTriangular)
                        @test rdiv!(copy(A1), A2)::t1 ≈ A1 / A2 ≈ M1 / M2
                        @test ldiv!(A2, copy(A1))::t1 ≈ A2 \ A1 ≈ M2 \ M1
                    end
                    if (uplo1 != uplo2 && elty1 == elty2 != Int && t2 != UnitLowerTriangular && t2 != UnitUpperTriangular)
                        @test lmul!(adjoint(A1), copy(A2)) ≈ A1' * A2 ≈ M1' * M2
                        @test lmul!(transpose(A1), copy(A2)) ≈ transpose(A1) * A2 ≈ transpose(M1) * M2
                        @test ldiv!(adjoint(A1), copy(A2)) ≈ A1' \ A2 ≈ M1' \ M2
                        @test ldiv!(transpose(A1), copy(A2)) ≈ transpose(A1) \ A2 ≈ transpose(M1) \ M2
                    end
                    if (uplo1 != uplo2 && elty1 == elty2 != Int && t1 != UnitLowerTriangular && t1 != UnitUpperTriangular)
                        @test rmul!(copy(A1), adjoint(A2)) ≈ A1 * A2' ≈ M1 * M2'
                        @test rmul!(copy(A1), transpose(A2)) ≈ A1 * transpose(A2) ≈ M1 * transpose(M2)
                        @test rdiv!(copy(A1), adjoint(A2)) ≈ A1 / A2' ≈ M1 / M2'
                        @test rdiv!(copy(A1), transpose(A2)) ≈ A1 / transpose(A2) ≈ M1 / transpose(M2)
                    end
                end
            end

            for eltyB in (Float32, Float64, BigFloat, ComplexF32, ComplexF64, Complex{BigFloat})
                # Only test methods for the same element type and a single combination of mixed element types
                # to avoid too much compilation
                if !(elty1 == eltyB || elty1 ∈ (ComplexF32, Int) || eltyB ∈ (ComplexF32, Int))
                    continue
                end

                B = convert(Matrix{eltyB}, (elty1 <: Complex ? real(A1) : A1) * fill(1.0, n, n))

                Tri = Tridiagonal(rand(eltyB, n - 1), rand(eltyB, n), rand(eltyB, n - 1))
                C = Matrix{promote_type(elty1, eltyB)}(undef, n, n)
                mul!(C, Tri, A1)
                @test C ≈ Tri * M1
                Tri = Tridiagonal(rand(eltyB, n - 1), rand(eltyB, n), rand(eltyB, n - 1))
                mul!(C, A1, Tri)
                @test C ≈ M1 * Tri

                # Triangular-dense Matrix/vector multiplication
                @test A1 * B[:, 1] ≈ M1 * B[:, 1]
                @test A1 * B ≈ M1 * B
                @test transpose(A1) * B[:, 1] ≈ transpose(M1) * B[:, 1]
                @test A1'B[:, 1] ≈ M1'B[:, 1]
                @test transpose(A1) * B ≈ transpose(M1) * B
                @test A1'B ≈ M1'B
                @test A1 * transpose(B) ≈ M1 * transpose(B)
                @test adjoint(A1) * transpose(B) ≈ M1' * transpose(B)
                @test transpose(A1) * adjoint(B) ≈ transpose(M1) * adjoint(B)
                @test A1 * B' ≈ M1 * B'
                @test B * A1 ≈ B * M1
                @test transpose(B[:, 1]) * A1 ≈ transpose(B[:, 1]) * M1
                @test B[:, 1]'A1 ≈ B[:, 1]'M1
                @test transpose(B) * A1 ≈ transpose(B) * M1
                @test transpose(B) * adjoint(A1) ≈ transpose(B) * M1'
                @test adjoint(B) * transpose(A1) ≈ adjoint(B) * transpose(M1)
                @test B'A1 ≈ B'M1
                @test B * transpose(A1) ≈ B * transpose(M1)
                @test B * A1' ≈ B * M1'
                @test transpose(B[:, 1]) * transpose(A1) ≈ transpose(B[:, 1]) * transpose(M1)
                @test B[:, 1]'A1' ≈ B[:, 1]'M1'
                @test transpose(B) * transpose(A1) ≈ transpose(B) * transpose(M1)
                @test B'A1' ≈ B'M1'

                if eltyB == elty1
                    @test mul!(similar(B), A1, B) ≈ M1 * B
                    @test mul!(similar(B), A1, adjoint(B)) ≈ M1 * B'
                    @test mul!(similar(B), A1, transpose(B)) ≈ M1 * transpose(B)
                    @test mul!(similar(B), adjoint(A1), adjoint(B)) ≈ M1' * B'
                    @test mul!(similar(B), transpose(A1), transpose(B)) ≈ transpose(M1) * transpose(B)
                    @test mul!(similar(B), transpose(A1), adjoint(B)) ≈ transpose(M1) * B'
                    @test mul!(similar(B), adjoint(A1), transpose(B)) ≈ M1' * transpose(B)
                    @test mul!(similar(B), adjoint(A1), B) ≈ M1' * B
                    @test mul!(similar(B), transpose(A1), B) ≈ transpose(M1) * B
                    # test also vector methods
                    B1 = vec(B[1, :])
                    @test mul!(similar(B1), A1, B1) ≈ M1 * B1
                    @test mul!(similar(B1), adjoint(A1), B1) ≈ M1' * B1
                    @test mul!(similar(B1), transpose(A1), B1) ≈ transpose(M1) * B1
                end
                #error handling
                Ann, Bmm, bm = A1, Matrix{eltyB}(undef, n + 1, n + 1), Vector{eltyB}(undef, n + 1)
                @test_throws DimensionMismatch lmul!(Ann, bm)
                @test_throws DimensionMismatch rmul!(Bmm, Ann)
                @test_throws DimensionMismatch lmul!(transpose(Ann), bm)
                @test_throws DimensionMismatch lmul!(adjoint(Ann), bm)
                @test_throws DimensionMismatch rmul!(Bmm, adjoint(Ann))
                @test_throws DimensionMismatch rmul!(Bmm, transpose(Ann))

                # ... and division
                @test A1 \ B[:, 1] ≈ M1 \ B[:, 1]
                @test A1 \ B ≈ M1 \ B
                @test transpose(A1) \ B[:, 1] ≈ transpose(M1) \ B[:, 1]
                @test A1' \ B[:, 1] ≈ M1' \ B[:, 1]
                @test transpose(A1) \ B ≈ transpose(M1) \ B
                @test A1' \ B ≈ M1' \ B
                @test A1 \ transpose(B) ≈ M1 \ transpose(B)
                @test A1 \ B' ≈ M1 \ B'
                @test transpose(A1) \ transpose(B) ≈ transpose(M1) \ transpose(B)
                @test A1' \ B' ≈ M1' \ B'
                Ann, bm = A1, Vector{elty1}(undef, n + 1)
                @test_throws DimensionMismatch Ann \ bm
                @test_throws DimensionMismatch Ann' \ bm
                @test_throws DimensionMismatch transpose(Ann) \ bm
                if t1 == UpperTriangular || t1 == LowerTriangular
                    @test_throws SingularException ldiv!(t1(zeros(elty1, n, n)), fill(eltyB(1), n))
                    @test_throws SingularException ldiv!(t1(zeros(elty1, n, n)), fill(eltyB(1), n, 2))
                    @test_throws SingularException rdiv!(fill(eltyB(1), n, n), t1(zeros(elty1, n, n)))
                end
                @test B / A1 ≈ B / M1
                @test B / transpose(A1) ≈ B / transpose(M1)
                @test B / A1' ≈ B / M1'
                @test transpose(B) / A1 ≈ transpose(B) / M1
                @test B' / A1 ≈ B' / M1
                @test transpose(B) / transpose(A1) ≈ transpose(B) / transpose(M1)
                @test B' / A1' ≈ B' / M1'

                # Error bounds
                !(elty1 in (BigFloat, Complex{BigFloat})) && !(eltyB in (BigFloat, Complex{BigFloat})) && errorbounds(A1, A1 \ B, B)
            end
        end
    end
end

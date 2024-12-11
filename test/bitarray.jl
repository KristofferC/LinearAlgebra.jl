using LinearAlgebra, Test, Random

tc(r1::NTuple{N,Any}, r2::NTuple{N,Any}) where {N} = all(x->tc(x...), [zip(r1,r2)...])
tc(r1::BitArray{N}, r2::Union{BitArray{N},Array{Bool,N}}) where {N} = true
tc(r1::SubArray{Bool,N1,BitArray{N2}}, r2::SubArray{Bool,N1,<:Union{BitArray{N2},Array{Bool,N2}}}) where {N1,N2} = true
tc(r1::Transpose{Bool,BitVector}, r2::Union{Transpose{Bool,BitVector},Transpose{Bool,Vector{Bool}}}) = true
tc(r1::T, r2::T) where {T} = true
tc(r1,r2) = false

# vectors size
const v1 = 260
# matrices size
const n1, n2 = 17, 20
# arrays size
const s1, s2, s3, s4 = 5, 8, 3, 7

bitcheck(b::BitArray) = Test._check_bitarray_consistency(b)
bitcheck(x) = true

function check_bitop_call(ret_type, func, args...; kwargs...)
    r2 = func(map(x->(isa(x, BitArray) ? Array(x) : x), args)...; kwargs...)
    r1 = func(args...; kwargs...)
    ret_type ≢ nothing && (@test isa(r1, ret_type) || @show ret_type, typeof(r1))
    @test tc(r1, r2)
    @test isequal(r1, r2)
    @test bitcheck(r1)
end
macro check_bit_operation(ex, ret_type)
    @assert Meta.isexpr(ex, :call)
    Expr(:call, :check_bitop_call, esc(ret_type), map(esc, ex.args)...)
end
macro check_bit_operation(ex)
    @assert Meta.isexpr(ex, :call)
    Expr(:call, :check_bitop_call, nothing, map(esc, ex.args)...)
end


b1 = bitrand(v1)
b2 = bitrand(v1)
@check_bit_operation dot(b1, b2) Int

b1 = bitrand(n1, n2)
@test_throws ArgumentError tril(b1, -n1 - 2)
@test_throws ArgumentError tril(b1, n2)
@test_throws ArgumentError triu(b1, -n1)
@test_throws ArgumentError triu(b1, n2 + 2)
for k in (-n1 - 1):(n2 - 1)
    @check_bit_operation tril(b1, k) BitMatrix
end
for k in (-n1 + 1):(n2 + 1)
    @check_bit_operation triu(b1, k) BitMatrix
end

for sz = [(n1,n1), (n1,n2), (n2,n1)], (f,isf) = [(tril,istril), (triu,istriu)]
    _b1 = bitrand(sz...)
    @check_bit_operation isf(_b1) Bool
    _b1 = f(bitrand(sz...))
    @check_bit_operation isf(_b1) Bool
end

b1 = bitrand(n1,n1)
b1 .|= copy(b1')
@check_bit_operation issymmetric(b1) Bool
@check_bit_operation ishermitian(b1) Bool

b1 = bitrand(n1)
b2 = bitrand(n2)
@check_bit_operation kron(b1, b2) BitVector

b1 = bitrand(s1, s2)
b2 = bitrand(s3, s4)
@check_bit_operation kron(b1, b2) BitMatrix

b1 = bitrand(v1)
@check_bit_operation diff(b1) Vector{Int}

b1 = bitrand(n1, n2)
@check_bit_operation diff(b1, dims=1) Matrix{Int}
@check_bit_operation diff(b1, dims=2) Matrix{Int}

b1 = bitrand(n1, n1)
@test ((svdb1, svdb1A) = (svd(b1), svd(Array(b1)));
        svdb1.U == svdb1A.U && svdb1.S == svdb1A.S && svdb1.V == svdb1A.V)
@test ((qrb1, qrb1A) = (qr(b1), qr(Array(b1)));
        Matrix(qrb1.Q) == Matrix(qrb1A.Q) && qrb1.R == qrb1A.R)

b1 = bitrand(v1)
@check_bit_operation diagm(0 => b1) BitMatrix

b1 = bitrand(v1)
b2 = bitrand(v1)
@check_bit_operation diagm(-1 => b1, 1 => b2) BitMatrix

b1 = bitrand(n1, n1)
@check_bit_operation diag(b1)

# matrix methods for evalpoly(X, p) = ∑ₖ Xᵏ⁻¹ p[k]

# non-inplace fallback for evalpoly(X, p)
function _evalpoly(X::AbstractMatrix, p)
    Base.require_one_based_indexing(p)
    p0 = isempty(p) ? Base.reduce_empty_iter(+, p) : p[end]
    Xone = one(X)
    S = Base.promote_op(*, typeof(Xone), typeof(Xone))(Xone) * p0
    for i = length(p)-1:-1:1
        S = X * S + @inbounds(p[i] isa AbstractMatrix ? p[i] : p[i] * I)
    end
    return S
end

_scalarval(x::Number) = x
_scalarval(x::UniformScaling) = x.λ

"""
    evalpoly!(Y::AbstractMatrix, X::AbstractMatrix, p)

Evaluate the matrix polynomial ``Y = \\sum_k X^{k-1} p[k]``, storing the result
in-place in `Y`, for the coefficients `p[k]` (a vector or tuple).  The coefficients
can be scalars, matrices, or [`UniformScaling`](@ref).

Similar to `evalpoly`, but may be more efficient by working more in-place.  (Some
allocations may still be required, however.)
"""
function evalpoly!(Y::AbstractMatrix, X::AbstractMatrix, p::Union{AbstractVector,Tuple})
    @boundscheck axes(Y,1) == axes(Y,2) == axes(X,1) == axes(X,2)
    Base.require_one_based_indexing(p)

    N = length(p)
    pN = iszero(N) ? Base.reduce_empty_iter(+, p) : p[N]
    if pN isa AbstractMatrix
        Y .= pN
    elseif N > 1 && p[N-1] isa Union{Number,UniformScaling}
        # initialize Y to p[N-1] I + X p[N], in-place
        Y .= X .* _scalarval(pN)
        for i in axes(Y,1)
            @inbounds Y[i,i] += p[N-1] * I
        end
        N -= 1
    else
        # initialize Y to one(Y) * pN in-place
        for i in axes(Y,1)
            for j in axes(Y,2)
                @inbounds Y[i,j] = zero(Y[i,j])
            end
            @inbounds Y[i,i] += one(Y[i,i]) * pN
        end
    end
    if N > 1
        Z = similar(Y) # workspace for mul!
        for i = N-1:-1:1
            mul!(Z, X, Y)
            if p[i] isa AbstractMatrix
                Y .= p[i] .+ Z
            else
                # Y = p[i] * I + Z, in-place
                Y .= Z
                for j in axes(Y,1)
                    @inbounds Y[j,j] += p[i] * I
                end
            end
        end
    end
    return Y
end

# fallback cases: call out-of-place _evalpoly
Base.evalpoly(X::AbstractMatrix, p::Tuple) = _evalpoly(X, p)
Base.evalpoly(X::AbstractMatrix, ::Tuple{}) = zero(one(X)) # dimensionless zero, i.e. 0 * X^0
Base.evalpoly(X::AbstractMatrix, p::AbstractVector) = _evalpoly(X, p)

# optimized in-place cases, limited to types like homogeneous tuples with length > 1
# where we can reliably deduce the output type (= type of X * p[2]),
# and restricted to StridedMatrix (for now) so that we can be more confident that this is a performance win:
Base.evalpoly(X::StridedMatrix{<:Number}, p::Tuple{T, T, Vararg{T}}) where {T<:Union{Number, UniformScaling}} =
    evalpoly!(similar(X, Base.promote_op(*, eltype(X), typeof(_scalarval(p[2])))), X, p)
Base.evalpoly(X::StridedMatrix{<:Number}, p::Tuple{AbstractMatrix{T}, AbstractMatrix{T}, Vararg{AbstractMatrix{T}}}) where {T<:Number} =
    evalpoly!(similar(X, Base.promote_op(*, eltype(X), T)), X, p)
Base.evalpoly(X::StridedMatrix{<:Number}, p::AbstractVector{<:Union{Number, UniformScaling}}) =
    length(p) < 2 ? _evalpoly(X, p) : evalpoly!(similar(X, Base.promote_op(*, eltype(X), typeof(_scalarval(p[begin+1])))), X, p)
Base.evalpoly(X::StridedMatrix{<:Number}, p::AbstractVector{<:AbstractMatrix{<:Number}}) =
    length(p) < 2 ? _evalpoly(X, p) : evalpoly!(similar(X, Base.promote_op(*, eltype(X), eltype(p[begin+1]))), X, p)

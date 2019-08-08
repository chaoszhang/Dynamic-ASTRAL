struct ClusterHash
	r1::UInt64
	r2::UInt64
end

function Base.:+(x::ClusterHash, y::ClusterHash)::ClusterHash
	return ClusterHash(x.r1 + y.r1, x.r2 + y.r2)
end

function Base.:-(x::ClusterHash, y::ClusterHash)::ClusterHash
	return ClusterHash(x.r1 - y.r1, x.r2 - y.r2)
end

function Base.:<(x::ClusterHash, y::ClusterHash)::Bool
	return x.r1 < y.r1 || (x.r1 == y.r1 && x.r2 < y.r2)
end

function Base.hash(x::ClusterHash)::UInt64
	return x.r1
end


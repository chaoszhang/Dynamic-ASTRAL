mutable struct NodeMutableProperties{TNode, ScoreType<:Real}
	bestScore::ScoreType
	bestSmallChild::TNode
	
	function NodeMutableProperties{TNode, ScoreType}() where {TNode, ScoreType}
		return new{TNode, ScoreType}(typemin(ScoreType))
	end
	function NodeMutableProperties{TNode, ScoreType}(node::TNode) where {TNode, ScoreType}
		return new{TNode, ScoreType}(typemin(ScoreType), node)
	end
	function NodeMutableProperties{TNode, ScoreType}(bestScore::ScoreType, bestSmallChild::TNode) where {TNode, ScoreType}
		return new{TNode, ScoreType}(bestScore, bestSmallChild)
	end
end

struct Node{ScoreType<:Real}
	smallChildren::Vector{Node{ScoreType}}
	clusterHash::ClusterHash
	mutableProperties::NodeMutableProperties{Node{ScoreType}, ScoreType}
	#parents::Vector{Node, 1}
	
	function Node{ScoreType}(smallChildren::Vector{Node{ScoreType}}, clusterHash::ClusterHash, mutableProperties::NodeMutableProperties{Node{ScoreType}, ScoreType}) where ScoreType
		return new{ScoreType}(smallChildren, clusterHash, mutableProperties)
	end
	function Node{ScoreType}(clusterHash::ClusterHash) where ScoreType
		node = new{ScoreType}(Node{ScoreType}[], clusterHash, NodeMutableProperties{Node{ScoreType}, ScoreType}())
		node.mutableProperties.bestSmallChild = node
		return node
	end
end

function Base.:(==)(x::Node{ScoreType}, y::Node{ScoreType})::Bool where ScoreType<:Real
	return x.clusterHash == y.clusterHash
end

function Base.:<(x::Node{ScoreType}, y::Node{ScoreType})::Bool where ScoreType<:Real
	return x.clusterHash < y.clusterHash
end

function Base.hash(x::Node{ScoreType})::UInt64 where ScoreType<:Real
	return Base.hash(x.clusterHash)
end

abstract type TreeNode{ScoreType<:Real} end

mutable struct TreeLeafNode{ScoreType<:Real} <: TreeNode{ScoreType}
	node::Node{ScoreType}
	#score::ScoreType
	cluster::BitSet
	parent::TreeNode{ScoreType}
	
	function TreeLeafNode{ScoreType}(node::Node{ScoreType}, cluster::BitSet, parent::TreeNode{ScoreType}) where ScoreType
		return new{ScoreType}(node, cluster, parent)
	end
	function TreeLeafNode{ScoreType}() where ScoreType
		return new{ScoreType}()
	end
end

mutable struct TreeInternalNode{ScoreType<:Real} <: TreeNode{ScoreType}
	node::Node{ScoreType}
	#score::ScoreType
	cluster::BitSet
	parent::TreeNode{ScoreType}
	leftChild::TreeNode{ScoreType}
	rightChild::TreeNode{ScoreType}
	
	function TreeInternalNode{ScoreType}(node::Node{ScoreType}, cluster::BitSet, parent::TreeNode{ScoreType}, leftChild::TreeNode{ScoreType}, rightChild::TreeNode{ScoreType}) where ScoreType
		return new{ScoreType}(node, cluster, parent, leftChild, rightChild)
	end
	function TreeInternalNode{ScoreType}() where ScoreType
		return new{ScoreType}()
	end
end

mutable struct TreeRootNode{ScoreType<:Real} <: TreeNode{ScoreType}
	node::Node{ScoreType}
	#score::ScoreType
	cluster::BitSet
	leftChild::TreeNode{ScoreType}
	rightChild::TreeNode{ScoreType}
	
	function TreeRootNode{ScoreType}(node::Node{ScoreType}, cluster::BitSet, leftChild::TreeNode{ScoreType}, rightChild::TreeNode{ScoreType}) where ScoreType
		return new{ScoreType}(node, cluster, leftChild, rightChild)
	end
	function TreeRootNode{ScoreType}() where ScoreType
		return new{ScoreType}()
	end
end

mutable struct TreeMutableProperties{ScoreType<:Real}
	score::ScoreType
end

struct Tree{ScoreType<:Real, nLeaves, nInternalNodes}
	root::TreeRootNode{ScoreType}
	leaves::NTuple{nLeaves, TreeLeafNode{ScoreType}}
	internalNodes::NTuple{nInternalNodes, TreeInternalNode{ScoreType}}
	mutableProperties::TreeMutableProperties{ScoreType}
end

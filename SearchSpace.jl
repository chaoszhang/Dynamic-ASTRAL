module SearchSpace

export initialize, randomTree, promotable, marginalBenefitOfPromotion, promoteNode, score

include("ClusterHash.jl")
include("Node.jl")
include("TreeNode.jl")
include("TripartitionScoreTask.jl")

function getNode(cluster::ClusterHash)::Node{ScoreType}
	if !haskey(nodeDict, cluster)
		nodeDict[cluster] = Node{ScoreType}(cluster)
	end
	return nodeDict[cluster]
end

function Base.:+(x::Node, y::Node)::Node
	if x > y
		x, y = y, x
	end
	return getNode(x.clusterHash + y.clusterHash)
end

function unsafeGetTripartitionScore(h1::ClusterHash, h2::ClusterHash)::ScoreType
	if h1 > h2
		h1, h2 = h2, h1
	end
	return scoreDict[h1, h2]
end

function promotable(n::TreeNode{TypeOfScore})::Bool where TypeOfScore<:Real
	return typeof(n) != TreeRootNode{TypeOfScore} && typeof(n.parent) != TreeRootNode{TypeOfScore}
end

function computeBenefitOfPromotion(n::TreeNode{TypeOfScore})::Nothing where TypeOfScore<:Real
	p = n.parent
	b = (p.leftChild === n ? p.rightChild : p.leftChild)
	g = p.parent
	u = (g.leftChild === p ? g.rightChild : g.leftChild)
	addTripartitionScoreTask(u.node.clusterHash, b.node.clusterHash, u.cluster, b.cluster)
	addTripartitionScoreTask(n.node.clusterHash, u.node.clusterHash + b.node.clusterHash, n.cluster, union(u.cluster, b.cluster))
end

function getBenefitOfPromotion(n::TreeNode{TypeOfScore})::TypeOfScore where TypeOfScore<:Real
	p = n.parent
	b = (p.leftChild === n ? p.rightChild : p.leftChild)
	g = p.parent
	u = (g.leftChild === p ? g.rightChild : g.leftChild)
	return unsafeGetTripartitionScore(u.node.clusterHash, b.node.clusterHash) + unsafeGetTripartitionScore(u.node.clusterHash + b.node.clusterHash, n.node.clusterHash) - unsafeGetTripartitionScore(b.node.clusterHash, n.node.clusterHash) - unsafeGetTripartitionScore(p.node.clusterHash, u.node.clusterHash)
end

function marginalBenefitOfPromotion(nodes::Vector{TreeNode{TypeOfScore}})::Vector{TypeOfScore} where TypeOfScore<:Real
	for n in nodes
		computeBenefitOfPromotion(n)
	end
	computeTripartitionScore()
	return [getBenefitOfPromotion(n) for n in nodes]
end

function promoteNode(tree::Tree{TypeOfScore, nLeaves, nInternalNodes}, n::TreeNode{TypeOfScore})::Nothing where {TypeOfScore<:Real, nLeaves, nInternalNodes}
	p = n.parent
	b = (p.leftChild === n ? p.rightChild : p.leftChild)
	g = p.parent
	u = (g.leftChild === p ? g.rightChild : g.leftChild)
	tree.mutableProperties.score += unsafeGetTripartitionScore(u.node.clusterHash, b.node.clusterHash) + unsafeGetTripartitionScore(u.node.clusterHash + b.node.clusterHash, n.node.clusterHash) - unsafeGetTripartitionScore(b.node.clusterHash, n.node.clusterHash) - unsafeGetTripartitionScore(p.node.clusterHash, u.node.clusterHash)
	
	n.parent = g
	u.parent = p
	if p.leftChild === n
		p.leftChild = u
	else
		p.rightChild = u
	end
	if g.leftChild === u
		g.leftChild = n
	else
		g.rightChild = n
	end
	p.node = b.node + u.node
	p.cluster = union(b.cluster, u.cluster)
	return
end

function score(tree::Tree{TypeOfScore, nLeaves, nInternalNodes})::TypeOfScore where {TypeOfScore<:Real, nLeaves, nInternalNodes}
	return tree.mutableProperties.score
end

function initialize(TypeOfScore::DataType, nTaxa::Int, inputFile::String)::Nothing
	global ScoreType = TypeOfScore
	global nLeaves = nTaxa - 1
	global singletons = [Node{ScoreType}(ClusterHash(rand(UInt64), rand(UInt64))) for i in 1:nLeaves]
	global scoreDict = IdDict{Tuple{ClusterHash, ClusterHash}, ScoreType}()
	global nodeDict = IdDict{ClusterHash, Node{ScoreType}}()
	global fullCluster = BitSet(0:nLeaves)
	initializeTripartitionScoreTask(inputFile, 32, (nLeaves ÷ 64) + 1)
	return
end

function randomTree()
	leaves = [TreeLeafNode{ScoreType}() for i in 1:nLeaves]
	internalNodes = [TreeInternalNode{ScoreType}() for i in 1:nLeaves - 2]
	curNodes = Array{TreeNode{ScoreType}}(leaves)
	for i in 1:nLeaves
		leaves[i].node = singletons[i]
		leaves[i].cluster = BitSet([i])
	end
	for i in 1:nLeaves - 2
		j = mod(rand(Int), length(curNodes)) + 1
		n1, curNodes[j] = curNodes[j], curNodes[length(curNodes)]
		internalNodes[i].leftChild = n1
		n1.parent = internalNodes[i]
		pop!(curNodes)
		
		k = mod(rand(Int), length(curNodes)) + 1
		n2, curNodes[k] = curNodes[k], curNodes[length(curNodes)]
		internalNodes[i].rightChild = n2
		n2.parent = internalNodes[i]
		curNodes[length(curNodes)] = internalNodes[i]
		
		internalNodes[i].node = n1.node + n2.node
		internalNodes[i].cluster = union(n1.cluster, n2.cluster)
		addTripartitionScoreTask(n1.node.clusterHash, n2.node.clusterHash, n1.cluster, n2.cluster)
	end
	root = TreeRootNode{ScoreType}(curNodes[1].node + curNodes[2].node, union(curNodes[1].cluster, curNodes[2].cluster), curNodes[1], curNodes[2])
	curNodes[1].parent = root
	curNodes[2].parent = root
	addTripartitionScoreTask(curNodes[1].node.clusterHash, curNodes[2].node.clusterHash, curNodes[1].cluster, curNodes[2].cluster)
	computeTripartitionScore()
	score = sum([unsafeGetTripartitionScore(tn.leftChild.node.clusterHash, tn.rightChild.node.clusterHash) for tn in internalNodes]) + unsafeGetTripartitionScore(root.leftChild.node.clusterHash, root.rightChild.node.clusterHash)
	return Tree{ScoreType, nLeaves, nLeaves - 2}(root, NTuple{nLeaves, TreeLeafNode{ScoreType}}(leaves), NTuple{nLeaves - 2, TreeInternalNode{ScoreType}}(internalNodes), TreeMutableProperties{ScoreType}(score))
end

end

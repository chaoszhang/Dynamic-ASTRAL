include("SearchSpace.jl")
using .SearchSpace
initialize(Int64, 14, "examplePolytree.txt")
t = randomTree()

while true
	println(score(t))
	ns = vcat([n for n in t.leaves if promotable(n)], [n for n in t.internalNodes if promotable(n)])
	margins = marginalBenefitOfPromotion(ns)
	if (max(margins...) <= 0)
		break
	end
	promoteNode(t, ns[argmax(margins)])
end

function resum(n::Union{SearchSpace.TreeRootNode{ScoreType}, SearchSpace.TreeInternalNode{ScoreType}}) where ScoreType
	return SearchSpace.unsafeGetTripartitionScore(n.leftChild.node.clusterHash, n.rightChild.node.clusterHash) + resum(n.leftChild) + resum(n.rightChild)
end

function resum(n::SearchSpace.TreeLeafNode{ScoreType}) where ScoreType
	return 0
end

resum(t.root)
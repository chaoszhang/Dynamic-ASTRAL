struct TripartitionScoreTask
	h1::ClusterHash
	h2::ClusterHash
	c1::BitSet
	c2::BitSet
end

function initializeTripartitionScoreTask(inputFile::String, batchSize::Int, bitSetLength::Int)::Nothing
	ccall((:init, "astral"), Cvoid, (Cstring,), inputFile)
	global tripartitionBatchSize = batchSize
	global tripartitionScoreTasks = Vector{TripartitionScoreTask}()
	global clusterBitSetLength = bitSetLength
	return
end

function computeTripartitionScore()::Nothing
	if length(tripartitionScoreTasks) == 0
		return
	end
	result = zeros(ScoreType, tripartitionBatchSize)
	b1 = zeros(UInt64, tripartitionBatchSize * clusterBitSetLength)
	b2 = zeros(UInt64, tripartitionBatchSize * clusterBitSetLength)
	b3 = zeros(UInt64, tripartitionBatchSize * clusterBitSetLength)
	for i in 1:length(tripartitionScoreTasks)
		t = tripartitionScoreTasks[i]
		c3 = setdiff(fullCluster, t.c1, t.c2)
		b1[(i - 1) * clusterBitSetLength + t.c1.offset + 1:(i - 1) * clusterBitSetLength + t.c1.offset + length(t.c1.bits)] = t.c1.bits
		b2[(i - 1) * clusterBitSetLength + t.c2.offset + 1:(i - 1) * clusterBitSetLength + t.c2.offset + length(t.c2.bits)] = t.c2.bits
		b3[(i - 1) * clusterBitSetLength + c3.offset + 1:(i - 1) * clusterBitSetLength + c3.offset + length(c3.bits)] = c3.bits
	end
	ccall((:batchCompute, "astral"), Cvoid, (Cint, Ptr{Culonglong}, Ptr{Culonglong}, Ptr{Culonglong}, Ptr{Culonglong}), length(tripartitionScoreTasks), result, b1, b2, b3)
	for i in 1:length(tripartitionScoreTasks)
		t = tripartitionScoreTasks[i]
		scoreDict[t.h1, t.h2] = result[i]
	end
	empty!(tripartitionScoreTasks)
	return
end

function addTripartitionScoreTask(h1::ClusterHash, h2::ClusterHash, c1::BitSet, c2::BitSet)::Nothing
	if h1 > h2
		h1, h2, c1, c2 = h2, h1, c2, c1
	end
	if haskey(scoreDict, (h1, h2))
		return
	end
	scoreDict[h1, h2] = typemin(ScoreType)
	n1 = getNode(h1)
	n2 = getNode(h2)
	ns = getNode(h1 + h2)
	push!(ns.smallChildren, n1)
	push!(tripartitionScoreTasks, TripartitionScoreTask(h1, h2, c1, c2))
	if length(tripartitionScoreTasks) == tripartitionBatchSize
		computeTripartitionScore()
	end
	return
end

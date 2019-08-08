#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include "WeightCalculationInstructionGenerator.cpp"
#include "x86intrin.h"

#define NUM_BITS_IN_WORD 64
#define BATCH_SIZE 32 

using namespace std;

extern "C" {
void init(char*);

void batchCompute(int, uint64_t*, uint64_t*, uint64_t*, uint64_t*);
}

__attribute__((vector)) __attribute__ ((always_inline)) inline uint64_t G(uint16_t xa16, uint16_t xc16, uint16_t ya16, uint16_t yb16, uint16_t yc16, uint16_t zb16){
	uint32_t ybm = xa16 - (uint16_t) 1, zbm1 = xa16 + yc16 - (uint16_t) 2, zbm2 = xc16 + ya16 - (uint16_t) 2, xa = xa16, yc = yc16, xc = xc16, ya = ya16, yb = yb16, zb = zb16;
	uint64_t xayc = xa * yc, xaycm = yb * ybm + zb * zbm1, xcya = xc * ya, xcyam = zb * zbm2;
	return xayc * xaycm + xcya * xcyam;
}

__attribute__((vector)) __attribute__ ((always_inline)) inline uint64_t G(uint16_t a, uint16_t b, uint16_t c, uint16_t d, 
		uint16_t e, uint16_t f, uint16_t g, uint16_t h, uint16_t i){
	uint16_t xa = a, xb = b, xc = c, ya = d, yb = e, yc = f, za = g, zb = h, zc = i;
	return G(xa, xc, ya, yb, yc, zb) + G(xb, xa, yb, yc, ya, zc) + G(xc, xb, yc, ya, yb, za);
}

struct Polytree{
	int n = 0, listSize = 0, queueSize = 0, nWord = 0;
	int* __restrict__ queue = nullptr;
	
	~Polytree(){
		if (queue) delete queue;
	}
	
	void compute(uint64_t* __restrict__ result, const uint64_t* __restrict__ b) const{
		__attribute__((aligned(64))) uint16_t (* __restrict__ lst)[3 * BATCH_SIZE] = new __attribute__((aligned(64))) uint16_t[listSize][3 * BATCH_SIZE]{};
		__attribute__((aligned(64))) uint64_t weight[BATCH_SIZE] = {};
		
		for (int i = 0; i < n; i++){
			for (int j = 0; j < 3 * BATCH_SIZE; j++){
				lst[i][j] = (b[j * nWord + i / NUM_BITS_IN_WORD] & (1LL << (i % NUM_BITS_IN_WORD))) ? 1 : 0;
			}
		}
		
		for (int i = 0; i < queueSize; i++){
			int cmd = queue[i];
			if (cmd >= 0) {
				int y = cmd, x1 = queue[++i], x2 = queue[++i], c = 2 * queue[++i];
				#pragma ivdep
				#pragma simd vectorlength(32)
				#pragma vector always
				#pragma vector aligned
				for (int k = 0; k < BATCH_SIZE; k++){
					weight[k] += c * G(lst[x1][k], lst[x1][k + BATCH_SIZE], lst[x1][k + 2 * BATCH_SIZE], 
						lst[x2][k], lst[x2][k + BATCH_SIZE], lst[x2][k + 2 * BATCH_SIZE],
						lst[y][k], lst[y][k + BATCH_SIZE], lst[y][k + 2 * BATCH_SIZE]);
				}
			}
			else{
				int cmd2 = queue[++i];
				if (cmd2 < 0){
					int d = ~cmd, y = ~cmd2;
					__attribute__((aligned(64))) uint64_t tempWeight[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint16_t sx0[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint16_t sx1[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint16_t sx2[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint32_t sxy0[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint32_t sxy1[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint32_t sxy2[BATCH_SIZE] = {};
					for (int j = i + 1; j <= i + d; j++){
						#pragma ivdep
						#pragma simd vectorlength(32)
						#pragma vector always
						#pragma vector aligned
						for (int k = 0; k < BATCH_SIZE; k++){
							uint16_t q0 = lst[queue[j]][k], q1 = lst[queue[j]][k + BATCH_SIZE], q2 = lst[queue[j]][k + 2 * BATCH_SIZE];
							uint32_t q0_32 = q0, q1_32 = q1, q2_32 = q2;
							sx0[k] += q0;
							sx1[k] += q1;
							sx2[k] += q2;
							sxy0[k] += q1_32 * q2_32;
							sxy1[k] += q2_32 * q0_32;
							sxy2[k] += q0_32 * q1_32;
						}
					}
					__attribute__((aligned(64))) uint64_t my0[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint64_t my1[BATCH_SIZE] = {};
					__attribute__((aligned(64))) uint64_t my2[BATCH_SIZE] = {};
					for (int j = i + 1; j <= i + d; j++){
						#pragma ivdep
						#pragma simd vectorlength(32)
						#pragma vector always
						#pragma vector aligned
						for (int k = 0; k < BATCH_SIZE; k++){
							uint16_t q0 = lst[queue[j]][k], q1 = lst[queue[j]][k + BATCH_SIZE], q2 = lst[queue[j]][k + 2 * BATCH_SIZE];
							uint32_t q0_32 = q0, q1_32 = q1, q2_32 = q2;
							uint32_t q0c2 = q0_32 * q0_32 - q0_32, q1c2 = q1_32 * q1_32 - q1_32, q2c2 = q2_32 * q2_32 - q2_32;
							uint64_t sx0_q0 = sx0[k] - q0, sx1_q1 = sx1[k] - q1, sx2_q2 = sx2[k] - q2;
							
							my0[k] += (uint64_t) q1c2 * sx2_q2 + (uint64_t) q2c2 * sx1_q1;
							my1[k] += (uint64_t) q2c2 * sx0_q0 + (uint64_t) q0c2 * sx2_q2;
							my2[k] += (uint64_t) q0c2 * sx1_q1 + (uint64_t) q1c2 * sx0_q0;
						}
					}
					#pragma ivdep
					#pragma simd vectorlength(32)
					#pragma vector always
					#pragma vector aligned
					for (int k = 0; k < BATCH_SIZE; k++){
						tempWeight[k] = 2 * (lst[y][k] * my0[k] + lst[y][k + BATCH_SIZE] * my1[k] + lst[y][k + 2 * BATCH_SIZE] * my2[k]);
					}
					for (int j = i + 1; j <= i + d; j++){
						#pragma ivdep
						#pragma simd vectorlength(32)
						#pragma vector always
						#pragma vector aligned
						for (int k = 0; k < BATCH_SIZE; k++){
							uint16_t q0 = lst[queue[j]][k], q1 = lst[queue[j]][k + BATCH_SIZE], q2 = lst[queue[j]][k + 2 * BATCH_SIZE];
							uint32_t q0_32 = q0, q1_32 = q1, q2_32 = q2;
							uint32_t sx0_q0 = sx0[k] - q0, sx1_q1 = sx1[k] - q1, sx2_q2 = sx2[k] - q2;
							uint32_t q0q1 = q0_32 * q1_32, q1q2 = q1_32 * q2_32, q2q0 = q2_32 * q0_32;
							uint64_t q0c2 = q0_32 * q0_32 - q0_32, q1c2 = q1_32 * q1_32 - q1_32, q2c2 = q2_32 * q2_32 - q2_32;
							tempWeight[k] += (uint64_t) (2 * (sx1_q1 * sx2_q2) - sxy0[k] + q1q2) * q0c2 
								+ (uint64_t) (2 * (sx2_q2 * sx0_q0) - sxy1[k] + q2q0) * q1c2
								+ (uint64_t) (2 * (sx0_q0 * sx1_q1) - sxy2[k] + q0q1) * q2c2;
						}
					}
					i += d + 1;
					int c = queue[i];
					#pragma ivdep
					#pragma simd vectorlength(32)
					#pragma vector always
					#pragma vector aligned
					for (int k = 0; k < BATCH_SIZE; k++){
						weight[k] += c * tempWeight[k];
					}
				}
				else {
					int y = ~cmd, x1 = cmd2, x2 = queue[++i];
					#pragma ivdep
					#pragma simd vectorlength(32)
					#pragma vector always
					#pragma vector aligned
					for (int k = 0; k < 3 * BATCH_SIZE; k++){
						lst[y][k] = lst[x1][k] + lst[x2][k];
					}
				}
			}
		}
		delete lst;
		#pragma ivdep
		#pragma simd vectorlength(32)
		#pragma vector always
		#pragma vector aligned
		for (int k = 0; k < BATCH_SIZE; k++){
			result[k] = weight[k];
		}//cerr << "result = " << result[0] << endl;
	}
	
} pt;

void init(char* cfile){
	string file(cfile);
	ifstream fin(file);
	fin >> pt.n;
	pt.nWord = (pt.n + 63) / 64;
	int numCluster, numSum, numScore;
	vector<int> idremap;
	for (int i = 0; i < pt.n; i++){
		string name;
		fin >> name;
		idremap.push_back(i);
	}
	fin >> numCluster >> numSum;
	pt.listSize = numCluster;
	vector<int> sumTarget, scoreTarget, scoreCnt;
	vector<vector<int> > sumSources, scoreSources;
	for (int i = 0; i < numSum; i++){
		int y, x0, x1;
		fin >> y >> x0 >> x1;
		sumTarget.push_back(y);
		sumSources.push_back({x0, x1});
	}
	fin >> numScore;
	for (int i = 0; i < numScore; i++){
		int y, x0, x1, c;
		fin >> y >> x0 >> x1 >> c;
		scoreTarget.push_back(y);
		scoreSources.push_back({x0, x1});
		scoreCnt.push_back(c);
	}
	WeightCalculationInstructionGenerator ig(idremap, sumTarget, sumSources, scoreTarget, scoreSources, scoreCnt, pt.n);
	const vector<int>& inst = ig.getInstruction();
	pt.queueSize = inst.size();
	pt.queue = new int[inst.size()];
	memcpy(pt.queue, inst.data(), inst.size() * sizeof(int));
}

void batchCompute(int size, uint64_t* jres, uint64_t* jb1, uint64_t* jb2, uint64_t* jb3){
	uint64_t result[BATCH_SIZE] = {};
	uint64_t* b = new uint64_t[3 * BATCH_SIZE * pt.nWord]{};
	for (int i = 0; i < size; i += BATCH_SIZE){
		int batchSize = min(BATCH_SIZE, size - i);
		memcpy(b + 0 * BATCH_SIZE * pt.nWord, jb1 + i * pt.nWord, batchSize * pt.nWord * sizeof(uint64_t));
		memcpy(b + 1 * BATCH_SIZE * pt.nWord, jb2 + i * pt.nWord, batchSize * pt.nWord * sizeof(uint64_t));
		memcpy(b + 2 * BATCH_SIZE * pt.nWord, jb3 + i * pt.nWord, batchSize * pt.nWord * sizeof(uint64_t));
		pt.compute(result, b);
		for (int j = 0; j < batchSize; j++){
			jres[i + j] = result[j];
		}
	}
	delete b;
}


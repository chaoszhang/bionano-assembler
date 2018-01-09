#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <algorithm>
#include <string>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic> 


using namespace std::chrono;
high_resolution_clock::time_point t1;
//duration<double> t2(0);

using namespace std;

#define MAX_SEED_LENGTH 5
#define MAX_SEED_SEARCH_SPACE_LENGTH 10
#define SEED_INTERVAL_SIZE 64
#define CHROMOSONE_SPACER 4294967296LL
#define BIN_SIZE 20

//#define DEBUG_MODE

mutex mtxIn, mtxOut;
atomic<bool> noMoreInput(false);

namespace GP{
	const long long seedLenVar = 250;
	const float seedLenModifier = 0.85;
	//const long long seedControlFactor = 100;
	//const long long seedMaxSize = 500000;
	const long long seedLenBoundary[SEED_INTERVAL_SIZE] = {366,734,1103,1472,1842,2213,2585,2959,3333,3710,4087,4466,4847,5230,5615,6002,6391,6784,7179,7577,7978,8383,8792,9206,9624,10047,10475,10910,11350,11798,12254,12717,13190,13672,14164,14668,15185,15714,16258,16818,17395,17990,18606,19243,19904,20590,21304,22048,22825,23637,24488,25381,26318,27305,28344,29442,30602,31830,33131,34512,35980,37542,39205,500000};
	const long long mergeThreshold = 1000;
	const long long mergeSpacer = 500;
	//const int FNThreshold = 1;
	const int seedBinVolumeThreshold = 2;
	const int seedFPThreshold = 1;
	const int FPThreshold = 1;
	const int DPMaxSearchSpace = 5;
	const float alignmentSizeDifferenceMaxTolerance = 5000;
	const int seedMinIntervalCnts = 1;
	const float seedingPortion = 0.6;
	const float FPScorePenalty = 1.7;
	const float FNScorePenalty = 1;
	const float oneIntervalBonusFactor = 1.2;
	const float twoSitesStartingFraction = 0.25;
	
	/*
		for (int i = 1; i < SEED_INTERVAL_SIZE; i++){
			seedLenBoundary[i - 1] = (seedLenVar + seedControlFactor / sqrt(seedLenModifier)) * i - seedControlFactor
				+ seedControlFactor * pow(sqrt(seedLenModifier), -i);
		}
		seedLenBoundary[SEED_INTERVAL_SIZE - 1] = seedMaxSize;
	*/
	
	inline int seedBinLocator(long long len){
		if (len >= seedLenBoundary[SEED_INTERVAL_SIZE - 1]) return -1;
		int low = 0, high = SEED_INTERVAL_SIZE;
		while (low + 1 != high){
			int mid = (low + high) / 2;
			if (len < seedLenBoundary[mid - 1]) high = mid;
			else low = mid;
		}
		return low;
	}
	
	inline float scoringFunction(float refLength, float seqLength){
		// h1 = k/w; h2 = ((2-2a) - k)/w; tail = e^(-ak/w); 
		const float k = 0.05, a = 0.025, delta = 750;
		float w = 0, x = 0;
		if (0.85 * refLength - delta < 0) {
			w = refLength + delta;
			x = 0;
		}
		else {
			w = 0.15 * refLength + 2 * delta;
			x = 0.85 * refLength - delta;
		}
		if (seqLength > refLength + delta) return log(k / w) - k / a * (seqLength - refLength - delta) / w;
		else if (seqLength < 0.85 * refLength - delta) return log(k / w) - k / a * (refLength * 0.85 - delta - seqLength) / w;
		else {
			if (seqLength > refLength) return log(((refLength + delta - seqLength) / delta * 2 * (1 - a - k) + k) / w);
			else return log(((seqLength - x) / (refLength - x) * 2 * (1 - a - k) + k) / w);
		}
	}
	
	inline float scoringFunction(float seqLength){
		if (seqLength < 100) return -log(100) - log(log(1000));
		else if (seqLength > 1e5) return -log(1e5) - log(log(1000));
		else return -log(seqLength) - log(log(1000));
	}
}

struct Genome{
	struct Seed{
		long long len[MAX_SEED_SEARCH_SPACE_LENGTH];
		int size, posID;
		
		Seed(long long* p, int size, int posID): size(size), posID(posID){
			for (int i = 0; i < size; i++){
				len[i] = p[i + 1] - p[i];
			}
		}
		
		void format(){
			int j = 0;
			for (int i = 0; i < size; i++){
				if (len[i] == 0) continue;
				len[j] = len[i];
				j++;
			}
			size = j;
		}
		
		int hashValue(){
			int hash = 1;
			for (int i = 0; i < size; i++){
				if  (len[i] == 0) continue;
				int t = GP::seedBinLocator(len[i]);
				if (t == -1) return 0;
				hash = hash * SEED_INTERVAL_SIZE + t;
			}
			return hash;
		}
		
		void FNSearchSpace(unordered_set<int> &s1, unordered_set<int> &s2, unordered_set<int> &s3, unordered_set<int> &s4){
			searchSpace(s1);
			searchSpace(s2);
			mergeSearchSpace(s3, 0);
			mergeSearchSpace(s4, 0);
			for (int i = 1; i < size; i++){
				Seed sd = *this;
				sd.len[i - 1] = len[i - 1] + len[i];
				sd.len[i] = 0;
				sd.format();
				sd.searchSpace(s2);
				sd.mergeSearchSpace(s4, 0);
			}
		}
		
		void mergeSearchSpace(unordered_set<int> &s, int i){
			if (i + 2 >= size) {
				Seed sd = *this;
				sd.searchSpace(s);
				return;
			}
			mergeSearchSpace(s, i + 1);
			if (len[i + 1] >= GP::mergeThreshold) return;
			Seed sd = *this;
			sd.len[i] = 0;
			for (int dist1 = 0; dist1 < len[i + 1]; dist1 += GP::mergeSpacer){
				sd.len[i + 1] = len[i] + dist1;
				sd.len[i + 2] = len[i + 1] - dist1 + len[i + 2];
				sd.mergeSearchSpace(s, i + 1);
			}
			sd.len[i + 1] = len[i] + len[i + 1];
			sd.len[i + 2] = len[i + 2];
			sd.mergeSearchSpace(s, i + 1);
		}
		
		void searchSpace(unordered_set<int> &s){
			format();
			if (size > MAX_SEED_LENGTH) return;
			int hash = hashValue();
			if (hash != 0) s.insert(hash);
		}
		
		void FNAddSearchSpace(unordered_multimap<int, int> &hash, unordered_set<long long> &presence, unordered_map<int, int> &c1, unordered_map<int, int> &c2, unordered_map<int, int> &c3, unordered_map<int, int> &c4){
			addSearchSpace(hash, presence, c1);
			addSearchSpace(hash, presence, c2);
			mergeAddSearchSpace(hash, presence, c3, 0);
			mergeAddSearchSpace(hash, presence, c4, 0);
			for (int i = 1; i < size; i++){
				Seed sd = *this;
				sd.len[i - 1] = len[i - 1] + len[i];
				sd.len[i] = 0;
				sd.format();
				sd.addSearchSpace(hash, presence, c2);
				sd.mergeAddSearchSpace(hash, presence, c4, 0);
			}
		}
		
		void mergeAddSearchSpace(unordered_multimap<int, int> &hash, unordered_set<long long> &presence, unordered_map<int, int> &c, int i){
			if (i + 2 >= size) {
				Seed sd = *this;
				sd.addSearchSpace(hash, presence, c);
				return;
			}
			mergeAddSearchSpace(hash, presence, c, i + 1);
			if (len[i + 1] >= GP::mergeThreshold) return;
			Seed sd = *this;
			sd.len[i] = 0;
			for (int dist1 = 0; dist1 < len[i + 1]; dist1 += GP::mergeSpacer){
				sd.len[i + 1] = len[i] + dist1;
				sd.len[i + 2] = len[i + 1] - dist1 + len[i + 2];
				sd.mergeAddSearchSpace(hash, presence, c, i + 1);
			}
			sd.len[i + 1] = len[i] + len[i + 1];
			sd.len[i + 2] = len[i + 2];
			sd.mergeAddSearchSpace(hash, presence, c, i + 1);
		}
		
		void addSearchSpace(unordered_multimap<int, int> &hash, unordered_set<long long> &presence, unordered_map<int, int> &c){
			format();
			if (size > MAX_SEED_LENGTH) return;
			long long h = hashValue();
			if (h == 0) return;
			if (c[h] > GP::seedBinVolumeThreshold) return;
			if (presence.count((h << 32) + posID)) return;
			hash.emplace(h, posID);
			presence.insert((h << 32) + posID);
		}
	};
	
	struct AlignmentResult{
		int startID = 0, endID = 0, qStart = 0, qEnd = 0;
		float absoluteAlignmentScore = 0, normalizedAlignmentScore = 0;
		
		bool operator < (const AlignmentResult &o) const{
			return absoluteAlignmentScore > o.absoluteAlignmentScore;
		}
		
		#ifdef DEBUG_MODE
		vector<int> ref, query;
		#endif
	};
	
	vector<long long> pos;
	vector<int> genomeID, locID, loc;
	vector<string> genomeName;
	vector<bool> strandReversed;
	unordered_multimap<int, int> hash;
	
	Genome(vector<vector<long long> > &positions, vector<string> &names){
		// unit tested
		unordered_map<int, int> blockCounter1, blockCounter2, blockCounter3, blockCounter4;
		unordered_set<int> hashset1, hashset2, hashset3, hashset4;
		unordered_set<long long> presence;
		
		for (int i = 0; i < positions.size(); i++){
			genomeName.push_back(names[i]);
			for (int j = 0; j < positions[i].size(); j++){
				pos.push_back(CHROMOSONE_SPACER * (2 * i) + positions[i][j]);
				genomeID.push_back(i);
				locID.push_back(j);
				loc.push_back(positions[i][j]);
				strandReversed.push_back(false);
			}
			for (int j = positions[i].size() - 1; j >= 0; j--){
				pos.push_back(CHROMOSONE_SPACER * (2 * i + 1) + positions[i].back() - positions[i][j]);
				genomeID.push_back(i);
				locID.push_back(j);
				loc.push_back(positions[i][j]);
				strandReversed.push_back(true);
			}
		}
		
		for (int i = 0; i < pos.size(); i++){
			hashset1.clear();
			hashset2.clear();
			hashset3.clear();
			hashset4.clear();
			for (int j = 1; j <= MAX_SEED_SEARCH_SPACE_LENGTH && i + j < pos.size(); j++){
				Seed s(&pos[i], j, i);
				s.FNSearchSpace(hashset1, hashset2, hashset3, hashset4);
			}
			for (int h: hashset1) if (blockCounter1[h] <= GP::seedBinVolumeThreshold) blockCounter1[h]++;
			for (int h: hashset2) if (blockCounter2[h] <= GP::seedBinVolumeThreshold) blockCounter2[h]++;
			for (int h: hashset3) if (blockCounter3[h] <= GP::seedBinVolumeThreshold) blockCounter3[h]++;
			for (int h: hashset4) if (blockCounter4[h] <= GP::seedBinVolumeThreshold) blockCounter4[h]++;
		}
		
		for (int i = 0; i < pos.size(); i++){
			for (int j = 1; j <= MAX_SEED_SEARCH_SPACE_LENGTH && i + j < pos.size(); j++){
				Seed s(&pos[i], j, i);
				s.FNAddSearchSpace(hash, presence, blockCounter1, blockCounter2, blockCounter3, blockCounter4);
			}
		}
	}
	
	void seedMapper(unordered_set<int> &h, vector<long long> &p, int i, int len, int nSkip, int hvalue) const{
		auto range = hash.equal_range(hvalue);
    	for (auto it = range.first; it != range.second; ++it){
    		h.insert(it->second);
    	}
		if (len == MAX_SEED_LENGTH) return;
		for (int j = 1; j <= nSkip + 1 && i + j < p.size(); j++){
			int lower = GP::seedBinLocator((p[i + j] - p[i]) - GP::seedLenVar);
			if (lower == -1) continue;
			int upper = GP::seedBinLocator((p[i + j] - p[i]) / GP::seedLenModifier + GP::seedLenVar);
			if (upper == -1) upper = SEED_INTERVAL_SIZE - 1;
			for (int v = lower; v <= upper; v++) seedMapper(h, p, i + j, len + 1, nSkip - j + 1, hvalue * SEED_INTERVAL_SIZE + v);
		}
	}
	
	vector<int> seedMapping(vector<long long> &p) const{
		vector<int> result;
		unordered_set<int> hashset;
		for (int i = 0; i < p.size(); i++){
			hashset.clear();
			seedMapper(hashset, p, i, 0, GP::seedFPThreshold, 1);
			for (int h: hashset) result.push_back(h);
		}
		return result;
	}
	
	void frontAlignment(long long pStart, long long pEnd, int iLocator, const vector<float> &seq, vector<AlignmentResult> &results) const{
		int iStart = iLocator, iEnd = iLocator;
		while (pos[iStart] > pStart && iStart != 0) iStart--;
		while (pos[iEnd] < pEnd && iEnd != pos.size() - 1) iEnd++;
		int m = iEnd - iStart + 1, n = seq.size(), bBest = 0, eBest = 0;
		float sBest = 0;
		vector<float> s(n * m), ref;
		vector<int> b(n * m);
		#ifdef DEBUG_MODE
		vector<int> bt(n * m);
		int btBest = 0;
		#endif
		for (int i = iStart; i <= iEnd; i++) ref.push_back(pos[i] - pos[iStart]);
		for (int p = 0; p < n; p++){
			for (int i = 0; i < m; i++){
				b[p * m + i] = p * m + i;
				#ifdef DEBUG_MODE
				bt[p * m + i] = p * m + i;
				#endif
			}
		}
		for (int p = 1; p < n; p++){
			float fpPenalty = 0;
			for (int q = p - 1; q >= 0 && p - q - 1 <= GP::FPThreshold; q--, fpPenalty += GP::FPScorePenalty){
				float refScore = GP::scoringFunction(seq[p] - seq[q]) + fpPenalty;
				float minRefDiff = seq[p] - seq[q] - GP::alignmentSizeDifferenceMaxTolerance;
				float maxRefDiff = (seq[p] - seq[q]) / GP::seedLenModifier + GP::alignmentSizeDifferenceMaxTolerance;
				for (int i = 1; i < m; i++){
					float totalPenalty = refScore;
					for (int j = i - 1; j >= 0; j--, totalPenalty += GP::FNScorePenalty){
						if (ref[i] - ref[j] < minRefDiff) continue;
						if (ref[i] - ref[j] > maxRefDiff) break;
						float curS = (GP::scoringFunction(ref[i] - ref[j], seq[p] - seq[q]) - totalPenalty);
						if (s[p * m + i] < s[q * m + j] + curS) {
							s[p * m + i] = s[q * m + j] + curS;
							b[p * m + i] = b[q * m + j];
							#ifdef DEBUG_MODE
							bt[p * m + i] = q * m + j;
							#endif
						}
					}
				}
			}
			for (int i = 0; i < m; i++){
				if (s[p * m + i] > sBest){
					sBest = s[p * m + i];
					bBest = b[p * m + i];
					eBest = p * m + i;
					#ifdef DEBUG_MODE
					btBest = bt[p * m + i];
					#endif
				}
			}
			if (sBest > results[p].absoluteAlignmentScore){
				results[p].absoluteAlignmentScore = sBest;
				results[p].startID = iStart + (bBest % m);
				results[p].endID = iStart + (eBest % m);
				results[p].qStart = bBest / m;
				results[p].qEnd = eBest / m;
			}
		}
		/*
		#ifdef DEBUG_MODE
		vector<int> rRef, rQuery;
		rRef.push_back((eBest % m) + iStart);
		rQuery.push_back(eBest / m);
		for (int k = eBest; k != bt[k]; k = bt[k]){
			rRef.push_back((bt[k] % m) + iStart);
			rQuery.push_back(bt[k] / m);
		}
		for (int i = rRef.size() - 1; i >= 0; i--){
			result.ref.push_back(rRef[i]);
			result.query.push_back(rQuery[i]);
		}
		#endif
		*/
	}
	
	void backAlignment(long long pStart, long long pEnd, int iLocator, const vector<float> &seq, vector<AlignmentResult> &results) const{
		int iStart = iLocator, iEnd = iLocator;
		while (pos[iStart] > pStart && iStart != 0) iStart--;
		while (pos[iEnd] < pEnd && iEnd != pos.size() - 1) iEnd++;
		int m = iEnd - iStart + 1, n = seq.size(), bBest = 0, eBest = 0;
		float sBest = 0;
		vector<float> s(n * m), ref;
		vector<int> b(n * m);
		#ifdef DEBUG_MODE
		vector<int> bt(n * m);
		int btBest = 0;
		#endif
		for (int i = iEnd; i >= iStart; i--) ref.push_back(pos[iEnd] - pos[i]);
		for (int p = 0; p < n; p++){
			for (int i = 0; i < m; i++){
				b[p * m + i] = p * m + i;
				#ifdef DEBUG_MODE
				bt[p * m + i] = p * m + i;
				#endif
			}
		}
		for (int p = 1; p < n; p++){
			float fpPenalty = 0;
			for (int q = p - 1; q >= 0 && p - q - 1 <= GP::FPThreshold; q--, fpPenalty += GP::FPScorePenalty){
				float refScore = GP::scoringFunction(seq[p] - seq[q]) + fpPenalty;
				float minRefDiff = seq[p] - seq[q] - GP::alignmentSizeDifferenceMaxTolerance;
				float maxRefDiff = (seq[p] - seq[q]) / GP::seedLenModifier + GP::alignmentSizeDifferenceMaxTolerance;
				for (int i = 1; i < m; i++){
					float totalPenalty = refScore;
					for (int j = i - 1; j >= 0; j--, totalPenalty += GP::FNScorePenalty){
						if (ref[i] - ref[j] < minRefDiff) continue;
						if (ref[i] - ref[j] > maxRefDiff) break;
						float curS = (GP::scoringFunction(ref[i] - ref[j], seq[p] - seq[q]) - totalPenalty);
						if (s[p * m + i] < s[q * m + j] + curS) {
							s[p * m + i] = s[q * m + j] + curS;
							b[p * m + i] = b[q * m + j];
							#ifdef DEBUG_MODE
							bt[p * m + i] = q * m + j;
							#endif
						}
					}
				}
			}
			for (int i = 0; i < m; i++){
				if (s[p * m + i] > sBest){
					sBest = s[p * m + i];
					bBest = b[p * m + i];
					eBest = p * m + i;
					#ifdef DEBUG_MODE
					btBest = bt[p * m + i];
					#endif
				}
			}
			if (sBest > results[n - 1 - p].absoluteAlignmentScore){
				results[n - 1 - p].absoluteAlignmentScore = sBest;
				results[n - 1 - p].startID = iEnd - (eBest % m);
				results[n - 1 - p].endID = iEnd - (bBest % m);
				results[n - 1 - p].qStart = n - 1 - eBest / m;
				results[n - 1 - p].qEnd = n - 1 - bBest / m;
			}
		}
		/*
		#ifdef DEBUG_MODE
		vector<int> rRef, rQuery;
		rRef.push_back((eBest % m) + iStart);
		rQuery.push_back(eBest / m);
		for (int k = eBest; k != bt[k]; k = bt[k]){
			rRef.push_back((bt[k] % m) + iStart);
			rQuery.push_back(bt[k] / m);
		}
		for (int i = rRef.size() - 1; i >= 0; i--){
			result.ref.push_back(rRef[i]);
			result.query.push_back(rQuery[i]);
		}
		#endif
		*/
	}
	
	void filterAlignmentSpace(unordered_map<long long, int> &alignmentSpace, const vector<int> seeds, long long step) const{
		if (seeds.size() == 0) return;
		int maxCount = 0;
		unordered_map<long long, pair<int, int> > intervalCnt;
		for (int s: seeds){
			pair<int, int> &e0 = intervalCnt[pos[s] / step];
			e0.first++;
			e0.second = s;
			if (e0.first > maxCount) maxCount = e0.first;
			pair<int, int> &e1 = intervalCnt[pos[s] / step - 1];
			e1.first++;
			e1.second = s;
			if (e1.first > maxCount) maxCount = e1.first;
		}
		int bin[BIN_SIZE + 1] = {}, sum = 0, cutoff = BIN_SIZE;
		for (pair<const long long, pair<int, int> > &e: intervalCnt){
			bin[e.second.first * BIN_SIZE / maxCount]++;
		}
		for (int i = BIN_SIZE; i >= 0; i--){
			sum += bin[i];
			if (sum <= GP::DPMaxSearchSpace) cutoff = i;
		}
		for (pair<const long long, pair<int, int> > &e: intervalCnt){
			if (e.second.first * BIN_SIZE / maxCount >= cutoff) {
				alignmentSpace[e.first] = e.second.second;
			}
		}
	}
	
	pair<AlignmentResult, AlignmentResult> align(const vector<long long> &seq) const{
		pair<AlignmentResult, AlignmentResult> result;
		unordered_map<long long, int> alignmentSpaceFront, alignmentSpaceBack;
		int frontMax, backMin;
		vector<long long> seqFront, seqBack;
		vector<float> seqFloat, seqFloatReverse;
		long long step = 2 * seq.back() / GP::seedLenModifier;
		for (frontMax = min((int) seq.size(), GP::seedMinIntervalCnts + 1); frontMax < seq.size(); frontMax++){
			if (seq[frontMax] > seq.back() * GP::seedingPortion) break;
		}
		for (int i = 0; i < frontMax; i++) seqFront.push_back(seq[i]);
		for (backMin = max(0, (int) seq.size() - GP::seedMinIntervalCnts - 1); backMin >= 0; backMin--){
			if (seq[backMin] < seq.back() * (1 - GP::seedingPortion)) break;
		}
		for (int i = backMin + 1; i < seq.size(); i++) seqBack.push_back(seq[i]);
		for (long long p: seq) seqFloat.push_back(p);
		for (int i = seq.size() - 1; i >= 0; i--) seqFloatReverse.push_back(seq[seq.size() - 1] - seq[i]);
		filterAlignmentSpace(alignmentSpaceFront, seedMapping(seqFront), step);
		filterAlignmentSpace(alignmentSpaceBack, seedMapping(seqBack), step);
		
		vector<AlignmentResult> front(seq.size()), back(seq.size());
		for (pair<const long long, int> &e: alignmentSpaceFront){
			frontAlignment(e.first * step, (e.first + 2) * step, e.second, seqFloat, front);
		}
		for (pair<const long long, int> &e: alignmentSpaceBack){
			backAlignment(e.first * step, (e.first + 2) * step, e.second, seqFloatReverse, back);
		}
		for (int i = 1; i < seq.size(); i++){
			double p = (i / (seq.size() - 1.0) + seq[i] / seq.back()) * 0.5;
			if (p < GP::twoSitesStartingFraction) continue;
			if (p > GP::twoSitesStartingFraction) break;
			if (front[i - 1].absoluteAlignmentScore + back[i].absoluteAlignmentScore
					> result.first.absoluteAlignmentScore + result.second.absoluteAlignmentScore){
				result.first = front[i - 1];
				result.second = back[i];
			}
		}
		if (front[seq.size() - 1].absoluteAlignmentScore > back[0].absoluteAlignmentScore){
			if (front[seq.size() - 1].absoluteAlignmentScore * GP::oneIntervalBonusFactor
					> result.first.absoluteAlignmentScore + result.second.absoluteAlignmentScore){
				result.first = front[seq.size() - 1];
				result.second = front[seq.size() - 1];
			}
		}
		else {
			if (back[0].absoluteAlignmentScore * GP::oneIntervalBonusFactor
					> result.first.absoluteAlignmentScore + result.second.absoluteAlignmentScore){
				result.first = back[0];
				result.second = back[0];
			}
		}
		return result;
	}
};

void job(const Genome &g, ifstream &fin, ofstream &fout){
	const int batchSize = 500;
	vector<vector<long long> > seqs;
	vector<long long> seq;
	int len;
	double p;
	string s;
	while (true){
		stringstream ss;
		mtxIn.lock();
		if (noMoreInput){
			mtxIn.unlock();
			return;
		}
		while (seqs.size() < batchSize){
			if (fin >> len){
				for (int i = 0; i < len; i++){
					fin >> p;
					seq.push_back(p);
				}
				if (len > 15 && len <= 40) seqs.push_back(seq);
				seq.clear();
			}
			else {
				noMoreInput = true;
				break;
			}
		}
		mtxIn.unlock();
		for (const vector<long long> &sq: seqs){
			auto r = g.align(sq);
			ss << sq.size() << "\t" << sq.back() - sq.front() << "\t"
				<< r.first.absoluteAlignmentScore << "\t" << r.second.absoluteAlignmentScore << "\t"
				<< g.genomeName[g.genomeID[r.first.startID]] << "\t" << g.loc[r.first.startID] << "\t" << g.loc[r.first.endID] << "\t"
				<< g.genomeName[g.genomeID[r.second.startID]] << "\t" << g.loc[r.second.startID] << "\t" << g.loc[r.second.endID] << "\t"
				<< sq[r.first.qStart] << "\t" << sq[r.first.qEnd] << "\t" << sq[r.second.qStart] << "\t" << sq[r.second.qEnd] << "\n";
		}
		seqs.clear();
		mtxOut.lock();
		while (getline(ss, s)) fout << s << endl;
		cerr << "time: " << duration_cast<duration<double> >(high_resolution_clock::now() - t1).count() << endl;
		mtxOut.unlock();
	}
}

int main(int argc, char** argv) {
	const int nThrd = 8;
	vector<thread> thrds;
	
	ifstream fin("data.txt");
	ifstream fref("ref.txt");
	ofstream fout("result.txt");
	
	vector<vector<long long> > pos;
	vector<string> names;
	string name;
	while (fref >> name){
		int n;
		long long loc;
		fref >> n;
		names.push_back(name);
		pos.emplace_back();
		for (int i = 0; i < n; i++){
			fref >> loc;
			pos.back().push_back(loc);
		}
	}
	
	const Genome g(pos, names);
	
	cerr << "Genome construction complete.\n";

	t1 = chrono::high_resolution_clock::now();
	for (int i = 1; i < nThrd; i++) thrds.emplace_back(job, ref(g), ref(fin), ref(fout));
	job(g, fin, fout);
	for (auto &e: thrds) e.join();
	return 0;
}

#include <bits/stdc++.h>
#include <random>

using namespace std;

// conditions
const int MAX_V = 500 + 10;
const int MAX_E = 5000 + 10;
const double TIME_LIMIT_SEC = 29.5 * 6;
const double FITNESS_CONSTANT = 5.5;

// parameters
const int NUMBER_OF_CHRS = 1024;
const int NUMBER_OF_NEW_CHRS = 128;
//const double GENERATION_GAP = 256. / NUMBER_OF_CHRS;

// Constants
const int PRIME = 13;
const int MOD = 1e9 + 7;

using ti = tuple<int, int, int>;

// #define DEBUG
// #define PRINT_DETAIL

int V, E;
int Edges[MAX_E][3];

int ChrIx = 0;
clock_t startsAt;
mt19937 gen;

// [0.0, 1.0)
double doubleRand() {
	return rand() * 1. / (RAND_MAX * 1. + 1);
}

// chromosome
struct CHR {
	vector<bool> genes;
	CHR() {}
	CHR(vector<bool> __genes): genes(__genes) {}
	void print(FILE *out) const {
		for(bool x: genes) fprintf(out, "%d", x);
		fprintf(out, "\n");
	}

	CHR crossover(const CHR &other) {
		// 1 ~ V-1
		int left = gen() % (V-2) + 1;
		int right = gen() % (V-2) + 1;
		if (right < left) swap(left, right);
		vector<bool> newGenes(V);
		copy(genes.begin(), genes.begin()+left, newGenes.begin());
		copy(other.genes.begin()+left, other.genes.begin()+right, newGenes.begin()+left);
		copy(genes.begin()+right, genes.end(), newGenes.begin()+right);
		CHR res = CHR(newGenes);
		return res;
	}
	void mutation(double probability) {
		if(probability <= 0) return;
		for(int i=0; i<V; i++) {
			if (doubleRand() < probability) {
				genes[i] = !genes[i];
			}
		}
	}
	int hash() const {
		long long res = 1;
		for(bool x: genes) {
			res = res * PRIME + x;
			res %= MOD;
		}
		return (int)res;
	}
};

int getCutValue(const CHR &chr) {
	int score = 0;
	for(int i=0; i<E; i++) {
		int a, b, w;
		a = Edges[i][0];
		b = Edges[i][1];
		w = Edges[i][2];
		if(chr.genes[a] != chr.genes[b]) {
			score += w;
		}
	}
	return score;
}

void getInput() {
	srand(time(NULL));
	startsAt = clock();
	random_device rd;
	gen = mt19937(rd());

	scanf("%d%d", &V, &E);
	vector<ti> temp;
	for(int i=0; i<E; i++) {
		int a, b, w;
		scanf("%d%d%d", &a, &b, &w);
		temp.emplace_back(a, b, w);
	}
	sort(temp.begin(), temp.end());
	for(int i=0; i<E; i++) {
		int a, b, w; tie(a, b, w) = temp[i];
		Edges[i][0] = a;
		Edges[i][1] = b;
		Edges[i][2] = w;
	}
}

// indexed chromosome
struct ICHR {
	int index;
	int score;
	int hash;
	CHR chr;

	ICHR(const CHR &chr, bool skip = false) {
		this->index = ChrIx++;
		if (not skip) {
			this->score = getCutValue(chr);
			this->hash = chr.hash();
		}
		this->chr = chr;
	}

	int getCost() const {
		return -score;
	}

	// minCost means bestCost, maxCost means worstCost
	double getFitness(int minCost, int maxCost) const {
		int cost = getCost();
		return 1. * (maxCost - cost) + 1. * (maxCost - minCost) / (FITNESS_CONSTANT - 1);
	}
};

CHR mkRandomChr() {
	vector<bool> res(V);
	for(int i=0; i<V; i++) {
		res[i] = rand() % 2;
	}
	return CHR(res);
}

bool operator<(const ICHR &a, const ICHR &b) {
	if (a.score != b.score) return a.score < b.score;
	// if (a.hash != b.hash) 
	return a.hash < b.hash;
	// return a.index < b.index;
}

// generation
struct GEN {
	int genCnt;
	set<ICHR> chrs;
	int sumOfCost;
	int numberOfChrs;

	GEN(int __numberOfChrs) {
		this->numberOfChrs = 0;
		genCnt = 0;
		chrs.clear();
		sumOfCost = 0;
		for(int i=0; i<__numberOfChrs; i++) {
			insertNewChr(ICHR(mkRandomChr()));
		}
	}

	ICHR randomSelection() {
		// just random O(n)
		int selectedIndex = gen() % numberOfChrs;
		auto it = chrs.begin();
		for(int i=0; i<selectedIndex; i++) {
			it++;
		}
		ICHR selected = *it;
		return selected;
	}

	ICHR fastRandomSelection() {
		ICHR random = ICHR(CHR(), true);
		auto it = chrs.begin();
		do {
			random.hash = gen() % MOD;
			int minCost = chrs.rbegin()->getCost(), maxCost = chrs.begin()->getCost();
			int cost = minCost + gen() % (maxCost - minCost + 1);
			random.score = -cost;
			//if (gen() % 1000 == 0) printf("score hash [%d %d] %d %d\n", -maxCost, -minCost, random.score, random.hash);
			it = chrs.lower_bound(random);
		}while(it == chrs.end());

		return *it;
	}

	ICHR selection() {
		// slow O(n)
		int minCost = chrs.rbegin()->getCost(), maxCost = chrs.begin()->getCost();
		double sumOfFitness = 1. * (maxCost - minCost) * numberOfChrs / (FITNESS_CONSTANT - 1) + maxCost * numberOfChrs - sumOfCost;
		if (minCost == maxCost) {
			return randomSelection();
		}

		double point = doubleRand() * sumOfFitness;
		double currentSum = 0;
		for(auto &chr: chrs) {
			currentSum += chr.getFitness(minCost, maxCost);
			if (currentSum >= point) {
				return chr;
			}
		}

		return *chrs.rbegin();
	}

	void replacementOne(ICHR &newChr) {
		eraseWorst();
		//eraseRandom();
		insertNewChr(newChr);
	}
	void eraseWorst() {
		auto it = chrs.begin();
		sumOfCost -= it->getCost();
		chrs.erase(it);
		numberOfChrs--;
	}
	void eraseRandom() {
		int selectedIndex = gen() % numberOfChrs;
		auto it = chrs.begin();
		for(int i=0; i<selectedIndex; i++) {
			it++;
		}
		sumOfCost -= it->getCost();
		chrs.erase(it);
		numberOfChrs--;
	}
	void insertNewChr(const ICHR &newChr) {
		chrs.insert(newChr);
		sumOfCost += newChr.getCost();
		numberOfChrs++;
	}

	bool isExistSameScoreAndHash(const ICHR &newChr) {
		auto it = chrs.lower_bound(newChr);
		if (it == chrs.end()) {
			return false;
		}
		if (it->score == newChr.score && it->hash == newChr.hash) {
			return true;
		}
		return false;
	}

	ICHR getBestIChr() {
		return *chrs.rbegin();
	}

	ICHR getWorstIChr() {
		return *chrs.begin();
	}

	double getAverageScore() {
		double result = 0;
		for(auto ichr: chrs) {
			result += ichr.score;
		}
		return result / numberOfChrs;
	}

	void print(FILE *out) {
		fprintf(out, "Generation %d:\n", genCnt);
		for(auto ichr: chrs) {
			fprintf(out, "[%d, %d]: ", ichr.score, ichr.index);
			ichr.chr.print(out);
		}
		int minCost, maxCost;
		double sumOfFitness = 0;
		minCost = maxCost = chrs.begin()->getCost();
		fprintf(out, "Costs ");
		for(auto chr: chrs) {
			minCost = min(minCost, chr.getCost());
			maxCost = max(maxCost, chr.getCost());
			fprintf(out, "%d ", chr.getCost());
		}
		fprintf(out, "\n");
		fprintf(out, "Fitnesses ");
		for(auto chr: chrs) {
			sumOfFitness += chr.getFitness(minCost, maxCost);
		}
		for(auto chr: chrs) {
			fprintf(out, "%f ", chr.getFitness(minCost, maxCost) / sumOfFitness);
		}
		fprintf(out, "\n");
	}
};

GEN Gen = GEN(0);

void process() {
	Gen = GEN(NUMBER_OF_CHRS);
	bool timeLimitExceeded = false;
	bool stopCondition = false;
	int maxScore = -MOD;
	int maxScoreNoProgressGenCount = 0;

	fprintf(stderr, "Gen\tHightest Function Value\tAverage Function Value\tLowest Function Value\tTime (s)");
#ifdef PRINT_DETAIL
	for(int i=0; i<NUMBER_OF_CHRS; i++) {
		fprintf(stderr, "\tchr%02d", i);
	}
#endif
	fprintf(stderr, "\n");
	do {
		clock_t current = clock();
		double duration = (current - startsAt) * 1. / CLOCKS_PER_SEC;

		int numberOfNews = NUMBER_OF_NEW_CHRS;
		vector<ICHR> newChrs;
		#ifdef DEBUG
		Gen.print(stderr);
		#endif
		for(int i=0; i<numberOfNews / 3 * 1; i++) {
			ICHR a = Gen.selection();
			ICHR b = Gen.selection();
			CHR c = a.chr.crossover(b.chr);
			ICHR d = ICHR(c);
			newChrs.push_back(d);
		}
		for(int i=0; i<numberOfNews / 3 * 2; i++) {
			CHR c = Gen.selection().chr;
			c.mutation(1. / V);
			ICHR d = ICHR(c);
			newChrs.push_back(d);
		}
		
		int numberOfInserted = 0;
		for(ICHR &newChr: newChrs) {
			if (Gen.isExistSameScoreAndHash(newChr) == false) {
				Gen.insertNewChr(newChr);
				numberOfInserted++;
			}
		}
		//*
		if (maxScoreNoProgressGenCount / 5 < Gen.numberOfChrs - NUMBER_OF_CHRS + numberOfInserted ||  Gen.numberOfChrs > NUMBER_OF_CHRS * 2) {
			for(int i=0; i<numberOfInserted; i++) {
				Gen.eraseWorst();
			}
		}
		//*/
		/*
		for(int i=0; i<numberOfInserted; i++) {
			Gen.eraseRandom();
		}
		//*/

		current = clock();
		duration = (current - startsAt) * 1. / CLOCKS_PER_SEC;
		timeLimitExceeded = duration > TIME_LIMIT_SEC;
		stopCondition = timeLimitExceeded;

		if (maxScore < Gen.getBestIChr().score) {
			maxScore = Gen.getBestIChr().score;
			maxScoreNoProgressGenCount = 0;
		}
		else {
			maxScoreNoProgressGenCount++;
		}

		if (Gen.genCnt % 50 == 0) {
			fprintf(stderr, "%d\t%d\t%f\t%d\t%f", Gen.genCnt, Gen.getBestIChr().score, -1., Gen.getWorstIChr().score, duration);
			//fprintf(stderr, "%d\t%d\t%f\t%d\t%f", Gen.genCnt, Gen.getBestIChr().score, Gen.getAverageScore(), Gen.getWorstIChr().score, duration);
			int cnt = 1;
			fprintf(stderr, "\n");
			for(auto it = Gen.chrs.rbegin(); it != Gen.chrs.rend(); it++) {
				auto ichr = *it;
				fprintf(stderr, "[%d %d] (%d)", ichr.score, ichr.hash, Gen.isExistSameScoreAndHash(ichr));
				cnt--;
				if (cnt == 0) break;
			}
			fprintf(stderr, "\n[Number: %d, Cnt: %d] (Inserted: %d)", Gen.numberOfChrs, maxScoreNoProgressGenCount, numberOfInserted);
			fprintf(stderr, "\n");
#ifdef PRINT_DETAIL
			for(auto &ichr: Gen.chrs) {
			//auto &ichr = *Gen.chrs.rbegin();
				fprintf(stderr, "\t");
				for(int i=0; i<N*2; i++) {
					int val = 0;
					for(int j=0; j<4; j++) {
						val *= 2;
						val += (int)(ichr.chr.genes[i*4+j]);
					}
					fprintf(stderr, "%X", val);
				}
			}
#endif
			fprintf(stderr, "\n");
		}
		Gen.genCnt++;
	}
	while(not stopCondition);
}

void printOutput() {
	ICHR best = Gen.getBestIChr();
	printf("%d ", best.score);
	for(int i=0; i<V; i++) {
		printf("%d", (int)best.chr.genes[i]);
	}
	puts("");
}

int main() {
	getInput();
	process();
	printOutput();
	return 0;
}

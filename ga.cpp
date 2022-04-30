#include <bits/stdc++.h>
#include <random>

using namespace std;

// conditions
const int MAX_V = 500 + 10;
const int MAX_E = 5000 + 10;
const double TIME_LIMIT_SEC = 29.5;
const double MUTATION_TIME_LIMIT_SEC = 28;
const double CHECK_STEADY_TIME_SEC = 29;
const double FITNESS_CONSTANT = 1.5;

// parameters
const int NUMBER_OF_CHRS = 256;
const double GENERATION_GAP = 4. / NUMBER_OF_CHRS;

using ti = tuple<int, int, int>;

// #define DEBUG
// #define PRINT_DETAIL

int V, E;
int Edges[MAX_E][3];

int ChrIx = 0;
clock_t startsAt;
mt19937 gen;
uniform_int_distribution<int> crossoverDistribution;
uniform_int_distribution<int> selectionDistribution = uniform_int_distribution<int>(0, NUMBER_OF_CHRS-1);

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
		int left = crossoverDistribution(gen);
		vector<bool> newGenes(V);
		copy(genes.begin(), genes.begin()+left, newGenes.begin());
		copy(other.genes.begin()+left, other.genes.end(), newGenes.begin()+left);
		CHR res = CHR(newGenes);
		#ifdef DEBUG
		printf("crossover %d\n", left);
		print(stderr);
		other.print(stderr);
		res.print(stderr);
		puts("");
		#endif
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

	crossoverDistribution = uniform_int_distribution<int>(1, V-1);
}

// indexed chromosome
struct ICHR {
	int index;
	int score;
	CHR chr;

	ICHR(const CHR &chr) {
		this->index = ChrIx++;
		this->score = getCutValue(chr);
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
	return a.index < b.index;
}

// generation
struct GEN {
	int genCnt;
	set<ICHR> chrs;
	int sumOfCost;

	GEN(int numberOfChrs) {
		genCnt = 0;
		chrs.clear();
		sumOfCost = 0;
		for(int i=0; i<numberOfChrs; i++) {
			insertNewChr(ICHR(mkRandomChr()));
		}
	}

	ICHR randomSelection() {
		// just random O(n)
		int selectedIndex = selectionDistribution(gen);
		auto it = chrs.begin();
		for(int i=0; i<selectedIndex; i++) {
			it++;
		}
		ICHR selected = *it;
		return selected;
	}

	ICHR selection() {
		// slow O(n)
		int minCost = chrs.rbegin()->getCost(), maxCost = chrs.begin()->getCost();
		double sumOfFitness = 1. * (maxCost - minCost) * NUMBER_OF_CHRS / (FITNESS_CONSTANT - 1) + maxCost * NUMBER_OF_CHRS - sumOfCost;
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
	}
	void eraseRandom() {
		int selectedIndex = selectionDistribution(gen);
		auto it = chrs.begin();
		for(int i=0; i<selectedIndex; i++) {
			it++;
		}
		chrs.erase(it);
	}
	void insertNewChr(const ICHR &newChr) {
		chrs.insert(newChr);
		sumOfCost += newChr.getCost();
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
		return result / NUMBER_OF_CHRS;
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
	bool mutationTimeLimitExceeded = false;
	bool steadyTimeLimitExceeded = false;
	bool allChromosomeSame = false;
	bool timeLimitExceeded = false;
	bool stopCondition = false;

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
		mutationTimeLimitExceeded = duration > MUTATION_TIME_LIMIT_SEC;

		int numberOfNews = GENERATION_GAP * NUMBER_OF_CHRS;
		vector<ICHR> newChrs;
		#ifdef DEBUG
		Gen.print(stderr);
		#endif
		for(int i=0; i<numberOfNews; i++) {
			ICHR a = Gen.selection();
			ICHR b = Gen.selection();
			CHR c = a.chr.crossover(b.chr);
			//double constant = 1;
			double constant = (TIME_LIMIT_SEC - duration) / 10;
            c.mutation(not mutationTimeLimitExceeded ? 1. / V  : 0);
			ICHR d = ICHR(c);
			newChrs.push_back(d);
		}
		for(ICHR &newChr: newChrs) {
			Gen.eraseWorst();
		}
		for(ICHR &newChr: newChrs) {
			Gen.insertNewChr(newChr);
		}

		current = clock();
		duration = (current - startsAt) * 1. / CLOCKS_PER_SEC;
		timeLimitExceeded = duration > TIME_LIMIT_SEC;
		steadyTimeLimitExceeded = duration > CHECK_STEADY_TIME_SEC;
		if (steadyTimeLimitExceeded) {
			bool isSame = true;
			for(auto &ichr: Gen.chrs) {
				if (ichr.chr.genes != Gen.chrs.begin()->chr.genes) {
					isSame = false;
					break;
				}
			}
			allChromosomeSame = true;
		}
		stopCondition = timeLimitExceeded or (steadyTimeLimitExceeded and allChromosomeSame);

		if (Gen.genCnt % 1000 == 0) {
			fprintf(stderr, "%d\t%d\t%f\t%d\t%f", Gen.genCnt, Gen.getBestIChr().score, Gen.getAverageScore(), Gen.getWorstIChr().score, duration);
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

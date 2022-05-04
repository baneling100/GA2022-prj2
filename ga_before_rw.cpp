#include <bits/stdc++.h>
#include <random>

using namespace std;

// conditions
const int MAX_V = 500 + 1;
const int MAX_E = 5000 + 1;
const double TIME_LIMIT_SEC = 29.5 * 6;
// const double FITNESS_CONSTANT = 5.5;

// parameters
const int INIT_NUMBER_OF_CHRS = 1024;
const int MAXIMUM_NUMBER_OF_CHRS = 1024 * 100;
const int NUMBER_OF_NEW_CHRS = 256;
//const double GENERATION_GAP = 256. / NUMBER_OF_CHRS;

// Constants
const int MOD = 1e9 + 7;

// #define DEBUG
// #define PRINT_DETAIL

int V, E;
int Edges[MAX_E][3];
vector<int> Ed[MAX_V];
int Renumber[MAX_V];
int RenumberInv[MAX_V];
int RenumberedEdges[MAX_E][3];

clock_t startsAt;
mt19937 gen;

int TwoPowerVMinusOne = 1;

inline int getComplementHash(int hash) {
	if (hash > TwoPowerVMinusOne)  return MOD - (hash - TwoPowerVMinusOne);
	return TwoPowerVMinusOne - hash;
}

// chromosome
class CHR {	
public:
	bool genes[MAX_V];
	CHR() {}
	CHR(bool __genes[MAX_V]) {
		for(int i=0; i<V; i++) {
			genes[i] = __genes[i];
		}
	}

	CHR* crossover(const CHR* other) {
		bool newGenes[MAX_V];
		// 0 ~ V-1
		int left = gen() % V;
		int right = gen() % V;
		int flip = gen() % 2;
		if (right < left) swap(left, right);
		for(int i=0; i<V; i++) {
			if (left <= i && i < right) newGenes[i] = flip ? !other->genes[i] : other->genes[i];
			else newGenes[i] = genes[i];
		}
		return new CHR(newGenes);
	}
	CHR* mutation() {
		CHR* newChr = new CHR(genes);
		int ix = gen() % V;
		newChr->genes[ix] = !newChr->genes[ix];
		return newChr;
	}
	int hash() const {
		int res = 0;
		for(int i=0; i<V; i++) {
			res = res * 2 + genes[i];
			if (res >= MOD) res -= MOD;
		}
		return min(res, getComplementHash(res));
	}
};

int getCutValue(const CHR* chr) {
	int score = 0;
	for(int i=0; i<E; i++) {
		int a, b, w;
		a = RenumberedEdges[i][0];
		b = RenumberedEdges[i][1];
		w = RenumberedEdges[i][2];
		if(chr->genes[a] != chr->genes[b]) {
			score += w;
		}
	}
	return score;
}

using ti = tuple<int, int, int>;

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
	// sort(temp.begin(), temp.end());
	for(int i=0; i<E; i++) {
		int a, b, w; tie(a, b, w) = temp[i];
		Edges[i][0] = a;
		Edges[i][1] = b;
		Edges[i][2] = w;
		Ed[a].push_back(b);
		Ed[b].push_back(a);
	}

	for(int i=0; i<V; i++) {
		TwoPowerVMinusOne *= 2;
		TwoPowerVMinusOne %= MOD;
	}
	TwoPowerVMinusOne += MOD - 1;
	TwoPowerVMinusOne %= MOD;
}

vector<int> degree(MAX_V, 0);
vector<bool> vis(MAX_V, false);
int renumberV = 0;
void dfs(int v, int p) {
	if(vis[v]) return;
	vis[v] = true;
	Renumber[v] = renumberV;
	RenumberInv[renumberV] = v;
	renumberV++;
	for(int w: Ed[v]) {
		dfs(w, v);
	}
}

using pi = pair<int, int>;
int st;
void renumber() {
	// TODO any more good algorithm ?
	int prenumber[MAX_V + 1];

	priority_queue<pi, vector<pi>, greater<pi>> Q;
	vector<int> visitOrder;
	st = gen() % V + 1;
	Q.emplace(-0, st);
	while(!Q.empty()) {
		int v = Q.top().second; Q.pop();
		if(vis[v]) continue;
		visitOrder.push_back(v);
		vis[v] = true;
		prenumber[v] = renumberV++;
		for(int w: Ed[v]) {
			if (vis[w]) continue;
			degree[w]++;
			Q.emplace(-degree[w], w);
		}
	}
	for(int i=1; i<=V; i++) {
		sort(Ed[i].begin(), Ed[i].end(), [&](int a, int b) {return prenumber[a] < prenumber[b];});
	}
	
	vis = vector<bool>(MAX_V, false);
	renumberV = 0;
	//for(int i: visitOrder) {
	for(int i=1; i<=V; i++) {
		dfs(i, -1);
	}

	// for(int i=1; i<=V; i++) {
	// 	Renumber[i] = i-1;
	// }

	for(int i=0; i<E; i++) {
		RenumberedEdges[i][0] = Renumber[Edges[i][0]];
		RenumberedEdges[i][1] = Renumber[Edges[i][1]];
		RenumberedEdges[i][2] = Edges[i][2];
	}

	// // TODO: DELETE
	// int test[MAX_V] = {-1, 0,72,40,191,33,226,219,120,232,172,63,158,221,216,280,195,130,296,128,189,32,112,8,89,249,148,30,48,207,31,162,87,186,205,55,250,255,56,291,271,53,90,214,251,147,51,113,104,233,151,47,201,223,260,54,256,157,150,194,6,70,242,9,39,265,295,266,111,139,110,171,294,187,142,64,263,135,241,102,28,282,62,292,199,239,170,15,125,132,276,167,23,253,269,49,138,203,259,52,163,97,240,179,287,197,145,185,192,180,208,21,161,212,181,283,272,60,154,206,85,183,285,45,65,160,143,248,68,37,80,2,289,215,140,4,159,231,24,224,127,178,177,156,38,274,227,277,184,261,136,267,281,213,58,106,222,286,95,220,7,41,73,275,247,230,173,293,22,204,190,46,254,35,200,19,126,16,108,202,155,188,175,114,252,270,149,81,235,105,5,118,119,264,218,82,198,124,176,168,29,100,96,44,210,123,117,14,268,115,238,278,133,36,61,121,134,17,20,237,57,169,146,288,79,3,71,77,234,246,13,94,245,25,290,66,217,273,228,193,78,164,116,92,257,59,43,34,75,103,284,50,211,26,144,153,243,129,244,141,101,12,93,182,69,262,88,76,236,258,131,109,229,1,152,137,165,166,209,86,279,18,225,107,27,98,99,174,122,91,74,10,84,42,11,83,196,67};
	// for(int i=0; i<E; i++) {
	// 	RenumberedEdges[i][0] = test[Edges[i][0]];
	// 	RenumberedEdges[i][1] = test[Edges[i][1]];
	// 	RenumberedEdges[i][2] = Edges[i][2];
	// }

	// TODO: DELETE
	// solution A / B
	// int testInv[MAX_V] = {4,6,7,13,14,20,29,34,40,43,52,59,67,73,94,97,108,110,113,119,131,133,151,153,156,159,169,170,179,181,185,194,204,208,236,239,252,278,1,2,3,5,8,9,10,11,12,15,16,17,18,19,21,22,23,24,25,26,27,28,30,31,32,33,35,36,37,38,39,41,42,44,45,46,47,48,49,50,51,53,54,55,56,57,58,60,61,62,63,64,65,66,68,69,70,71,72,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,95,96,98,99,100,101,102,103,104,105,106,107,109,111,112,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,132,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,152,154,155,157,158,160,161,162,163,164,165,166,167,168,171,172,173,174,175,176,177,178,180,182,183,184,186,187,188,189,190,191,192,193,195,196,197,198,199,200,201,202,203,205,206,207,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,237,238,240,241,242,243,244,245,246,247,248,249,250,251,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297};
	// int test[MAX_V];
	// for(int i=0; i<V; i++) test[testInv[i]] = i;
	// for(int i=0; i<E; i++) {
	// 	RenumberedEdges[i][0] = test[Edges[i][0]];
	// 	RenumberedEdges[i][1] = test[Edges[i][1]];
	// 	RenumberedEdges[i][2] = Edges[i][2];
	// }
}


CHR* mkRandomChr() {
	bool res[MAX_V];
	for(int i=0; i<V; i++) {
		res[i] = gen() % 2;
	}
	return new CHR(res);
}

// indexed chromosome
struct ICHR {
	int score;
	int hash;
	CHR* chr;

	ICHR() {
	}
	ICHR(CHR* chr) {
		this->score = getCutValue(chr);
		this->hash = chr->hash();
		this->chr = chr;
	}
};

bool operator<(const ICHR &a, const ICHR &b) {
	if (a.score != b.score) return a.score > b.score;
	return a.hash < b.hash;
}

// generation
struct GEN {
	int genCnt;
	int numberOfChrs;
	int bestScore;
	int bestScoreKeptGenCount;
	double bestScoreLastUpdatedAt = 0;
	CHR *chrs[MAXIMUM_NUMBER_OF_CHRS], *newChrs[NUMBER_OF_NEW_CHRS];
	// decresing order by score;
	ICHR ichrs[MAXIMUM_NUMBER_OF_CHRS + NUMBER_OF_NEW_CHRS];

	GEN(int initNumberOfChrs) {
		this->numberOfChrs = initNumberOfChrs;
		this->bestScore = -MOD;
		this->bestScoreKeptGenCount = 0;
		for(int i=0; i<numberOfChrs; i++) {
			chrs[i] = mkRandomChr();
			ichrs[i] = ICHR(chrs[i]);
		}
		sort(ichrs, ichrs + numberOfChrs);
		for(int i=0; i<numberOfChrs; i++) {
			chrs[i] = ichrs[i].chr;
		}
	}

	double averageScore() {
		double sum = 0;
		for(int i=0; i<numberOfChrs; i++) sum += ichrs[i].score;
		return sum / numberOfChrs;
	}

	void replace() {
		for(int i=0; i<NUMBER_OF_NEW_CHRS; i++) {
			ichrs[numberOfChrs + i] = ICHR(newChrs[i]);
		}
		int totalChrs = numberOfChrs + NUMBER_OF_NEW_CHRS;
		sort(ichrs, ichrs + totalChrs);

		int ichrsCnt = 1;
		for(int i=1; i<totalChrs; i++) {
			if (ichrs[i].score == ichrs[i-1].score && ichrs[i].hash == ichrs[i-1].hash)
				delete ichrs[i].chr;
			else
				ichrs[ichrsCnt++] = ichrs[i];
		}

		// if best score not incresing, increse number of chrs.
		if (ichrsCnt <= MAXIMUM_NUMBER_OF_CHRS && numberOfChrs - INIT_NUMBER_OF_CHRS < bestScoreKeptGenCount * 2) {
			numberOfChrs = ichrsCnt;
			for(int i=0; i<numberOfChrs; i++) {
				chrs[i] = ichrs[i].chr;
			}
			return;
		}

		// TODO just delete worst ?
		// int added = ichrsCnt - numberOfChrs;
		// for(int i=numberOfChrs - added; i<numberOfChrs; i++) {
		// 	swap(ichrs[i], ichrs[i+added]);
		// }
		
		// if it is too young, rescue it.
		// if (genCnt >= 1000) {
		// 	int swapIx = numberOfChrs - 1;
		// 	for(int i=ichrsCnt-1; i>=numberOfChrs; i--) {
		// 		if (ichrs[i].age >= genCnt - 100) {
		// 			swap(ichrs[swapIx--], ichrs[i]);
		// 		}
		// 	}
		// }

		for(int i=0; i<numberOfChrs; i++) {
			chrs[i] = ichrs[i].chr;
		}
		for(int i=numberOfChrs; i<ichrsCnt; i++) {
			delete ichrs[i].chr;
		}

		if (bestScore == ichrs[0].score)
			bestScoreKeptGenCount++;
		else {
			clock_t current = clock();
			bestScoreLastUpdatedAt = (current - startsAt) * 1. / CLOCKS_PER_SEC;
			bestScoreKeptGenCount = 0;
		}
		bestScore = max(bestScore, ichrs[0].score);
	}
};

GEN Gen = GEN(0);

void process() {
	Gen = GEN(INIT_NUMBER_OF_CHRS);
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
		int numberOfNews = NUMBER_OF_NEW_CHRS;
		#ifdef DEBUG
		Gen.print(stderr);
		#endif
		int crossoverNr = numberOfNews / 4 * 1;
		for(int i=0; i<crossoverNr; i++) {
			int a = gen() % Gen.numberOfChrs;
			int b = gen() % Gen.numberOfChrs;
			Gen.newChrs[i] = Gen.chrs[a]->crossover(Gen.chrs[b]);
		}
		for(int i=crossoverNr; i<numberOfNews; i++) {
			int a = gen() % Gen.numberOfChrs;
			Gen.newChrs[i] = Gen.chrs[a]->mutation();
		}

		Gen.replace();

		clock_t current = clock();
		double duration = (current - startsAt) * 1. / CLOCKS_PER_SEC;
		timeLimitExceeded = duration > TIME_LIMIT_SEC;
		stopCondition = timeLimitExceeded;

		if (Gen.genCnt % 50 == 0) {
			fprintf(stderr, "%d\t%d\t%f\t%d\t%f\t%d\n", Gen.genCnt, Gen.ichrs[0].score, Gen.averageScore(), Gen.ichrs[Gen.numberOfChrs-1].score, duration, Gen.numberOfChrs);
			for(int i=0; i<5; i++) {
				fprintf(stderr, "(%d %d) ", Gen.ichrs[i].score, Gen.ichrs[i].hash);
			}
			fprintf(stderr, "\n");
		}
		Gen.genCnt++;
	}
	while(not stopCondition);
}

void printOutput() {
	ICHR best = Gen.ichrs[0];

	// TODO: delete before submission
	int a = (int)(Gen.bestScoreLastUpdatedAt); double b = Gen.bestScoreLastUpdatedAt - a;
	printf("%4d [Gen: %7d St: %4d Last Updated At: %3d + %3f] ", best.score, Gen.genCnt, st, a, b);
	bool answer[MAX_V];
	for(int i=0; i<V; i++) {
		answer[RenumberInv[i]] = best.chr->genes[i];
	}
	for(int i=1; i<=V; i++) {
		printf("%d", (int)answer[i]);
	}
	printf(" ");

	bool isFirst = true;
	for(int i=1; i<=V; i++) {
		if (answer[i]) {
			if(not isFirst) printf(" ");
			if(isFirst) isFirst = false;
			printf("%d", i);
		}
	}

	puts("");
}

int main() {
	getInput();
	renumber();
	process();
	printOutput();
	return 0;
}

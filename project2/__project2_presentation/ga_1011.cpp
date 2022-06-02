#include <bits/stdc++.h>
#include <random>

using namespace std;

// conditions
const int MAX_V = 500 + 1;
const int MAX_E = 5000 + 1;
const double TIME_LIMIT_SEC = 29.5 * 6;
const double FITNESS_CONSTANT = 10;
const int FITNESS_GROUP = 8;

// parameters
const int INIT_NUMBER_OF_CHRS = 1024 * 10;
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

// [0, 1]
inline double getRandDouble() {
   return gen() * 1. / UINT_MAX;
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
   st = gen() % V + 1;
   Q.emplace(-0, st);
   while(!Q.empty()) {
      int v = Q.top().second; Q.pop();
      if(vis[v]) continue;
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
   for(int i=1; i<=V; i++) {
      dfs(i, -1);
   }
   for(int i=0; i<E; i++) {
      RenumberedEdges[i][0] = Renumber[Edges[i][0]];
      RenumberedEdges[i][1] = Renumber[Edges[i][1]];
      RenumberedEdges[i][2] = Edges[i][2];
   }
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
   double groupProbability[FITNESS_GROUP];

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
      double psum = 0;
      for(int i=0; i<FITNESS_GROUP; i++) {
         groupProbability[i] = 1 + (FITNESS_CONSTANT - 1) / (FITNESS_GROUP - 1) * (FITNESS_GROUP - 1 - i);
         psum += groupProbability[i];
      }
      for(int i=0; i<FITNESS_GROUP; i++) {
         groupProbability[i] /= psum;
         if (i >= 1) groupProbability[i] += groupProbability[i-1];
      }
   }

   double averageScore() {
      double sum = 0;
      for(int i=0; i<numberOfChrs; i++) sum += ichrs[i].score;
      return sum / numberOfChrs;
   }

   int selection() {
      double r = getRandDouble();
      int gix = FITNESS_GROUP - 1;
      for(int i=0; i<FITNESS_GROUP; i++) {
         if (r < groupProbability[i]) {
            gix = i;
            break;
         }
      }
      int minV = numberOfChrs / FITNESS_GROUP * gix;
      int maxV = min(numberOfChrs, numberOfChrs / FITNESS_GROUP * (gix+1));
      
      return minV + gen() % (maxV - minV);
   }

   void replace() {
      for(int i=0; i<NUMBER_OF_NEW_CHRS; i++) {
         ichrs[numberOfChrs + i] = ICHR(newChrs[i]);
      }
      int totalChrs = numberOfChrs + NUMBER_OF_NEW_CHRS;
      sort(ichrs, ichrs + totalChrs);

      int ichrsCnt = 1;
      for(int i=1; i<totalChrs; i++) {
         if (false && ichrs[i].score == ichrs[i-1].score && ichrs[i].hash == ichrs[i-1].hash)
            delete ichrs[i].chr;
         else
            ichrs[ichrsCnt++] = ichrs[i];
      }

      // // if best score not incresing, increse number of chrs.
      // if (ichrsCnt <= MAXIMUM_NUMBER_OF_CHRS && numberOfChrs - INIT_NUMBER_OF_CHRS < bestScoreKeptGenCount * 2) {
      //    numberOfChrs = ichrsCnt;
      //    for(int i=0; i<numberOfChrs; i++) {
      //       chrs[i] = ichrs[i].chr;
      //    }
      //    return;
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

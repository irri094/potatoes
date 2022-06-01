/** Mat Intersection per increment O(r*n) Weighted
Mat Intersection: per increment O(r^2*n); Evrythng
0-indexed erase base Mat to get better runtime*/
#include<bits/stdc++.h>
using namespace std;
typedef pair<int,int>PII; typedef vector<int>VI;
typedef vector<bool>VB;typedef long long CostType;
const CostType INF=1e18; const int BITSET_BITS=60;
struct Graph {
  vector<VI>edg;
  Graph(int nodes) : edg(nodes) { }
  void addEdge(int u,int v){ edg[u].push_back(v);}
  void clearGraph(){
    for (int i=0; i<edg.size();i++)edg[i].clear();
  }
};
struct Mat {
  virtual void updTknElmnt(const VB  &tkn) = 0;
  virtual bool canTkElmnt(int e) = 0;
  virtual bool canEx(int rmv, int ins) = 0;
};
struct ColorMat : Mat {
  VI elmntCol, canTakeAtMost, curTaking;
  int elmnts, clrs;
  ColorMat(int elmnts, int clrs):elmnts(elmnts),
    clrs(clrs), canTakeAtMost(clrs,1),
    curTaking(clrs,0) {
    elmntCol.reserve(elmnts);
  }
  void updTknElmnt(const VB &tkn) {
    fill(curTaking.begin(), curTaking.end(), 0);
    for (int i = 0; i < elmnts; i++) if (tkn[i]) {
      curTaking[ elmntCol[i] ]++;
    }
  }
  bool canTkElmnt(int e) {
    int col = elmntCol[e];
    return curTaking[col] != canTakeAtMost[col];
  }
  bool canEx(int rmv, int ins) {
    int colr = elmntCol[rmv],coli = elmntCol[ins];
    if (coli == colr) return true;
    return curTaking[coli] != canTakeAtMost[coli];
  }
};
struct GraphMat : Mat {
  vector< PII >ajs; int elmnts, grSz;
  GraphMat(int elmnts, int grSz) : forest(grSz) {
    this->elmnts = elmnts; this->grSz = grSz;
    ajs.reserve(elmnts);
  }
  Graph forest; VI start, finish, root;
  void treeDFS(int u, int p, int &tym) {
    start[u] = ++tym;
    for (int v : forest.edg[u]) {
      if (v == p) continue;
      root[v] = root[u]; treeDFS(v, u, tym);
    }
    finish[u] = tym;
  }
  bool inSubtree(int u, int x) {
  return start[u]<=start[x]&&finish[x]<=finish[u];
  }
  void updTknElmnt(const VB &tkn) {
    forest.clearGraph();
    for (int i = 0; i < elmnts; i++) if (tkn[i]) {
      forest.addEdge(ajs[i].first, ajs[i].second);
      forest.addEdge(ajs[i].second, ajs[i].first);
    }
    root=VI(grSz, -1); finish=start=VI(grSz,0);
    int tym = -1;
    for (int i = 0; i < grSz; i++)if(root[i]==-1){
      root[i] = i; treeDFS(i, -1, tym);
    }
  }
  bool canTkElmnt(int e) {
  return root[ajs[e].first]!=root[ajs[e].second];
  }
  bool canEx(int rmv, int ins) {
    if (canTkElmnt(ins)) return true;
    int u = ajs[rmv].first, p = ajs[rmv].second;
    if (start[p] > start[u]) swap(u, p);
    int x = ajs[ins].first, y = ajs[ins].second;
    return inSubtree(u, x) != inSubtree(u, y);
  }
};
struct BinMat : Mat {
  typedef bitset< BITSET_BITS > BitSet;
  struct Basis {
    int bitCnt;
    vector< BitSet >rdcd, combi;
    Basis(int bitCnt):bitCnt(bitCnt),rdcd(bitCnt),
      combi(bitCnt) {
      assert(BITSET_BITS == bitCnt);
    }
    void clearAll() {
      for (int i = 0; i < bitCnt; i++) {
        rdcd[i].reset(); combi[i].reset();
      }
    }
    BitSet canBeBuiltWith(BitSet x) {
      BitSet rt;
      for (int i = bitCnt-1; i >= 0; i--)
          if (x.test(i)) {
        if (!rdcd[i].test(i)) return BitSet();
        x ^= rdcd[i]; rt ^= combi[i];
      }
      return rt;
    }
    int addVector(BitSet x) {
      BitSet cm;
      for (int i = bitCnt-1; i >= 0; i--)
        if (x.test(i)) {
        if (!rdcd[i].test(i)) {
          rdcd[i]= x;combi[i]= cm.set(i);return i;
        } else {
          x ^= rdcd[i]; cm ^= combi[i];
        }
      }
      return -1;
    }
  };
  vector<BitSet> rows;int elmnts, bitCnt;
  BinMat(int elmnts, int bitCnt) : elmnts(elmnts),
    bitCnt(bitCnt), curBas(bitCnt), cycle(elmnts),
    rowMap(elmnts) { rows.reserve(elmnts);
  }
  vector< BitSet >cycle;
  VI rowMap;
  Basis curBas;
  void updTknElmnt(const VB &tkn) {
    curBas.clearAll();
    for (int i = 0; i < elmnts; i++) if (tkn[i]) {
      rowMap[i] = curBas.addVector(rows[i]);
    }
    for (int i = 0; i < elmnts; i++) if (!tkn[i]){
      cycle[i] = curBas.canBeBuiltWith(rows[i]);
    }
  }
  bool canTkElmnt(int e){ return !cycle[e].any();}
  bool canEx(int rmv, int ins) {
    if (canTkElmnt(ins)) return true;
    return cycle[ins].test( rowMap[rmv] );
  }
};
struct GrafDual : Mat {
  struct Bridge {
    vector< VI >adj; vector< PII >ajs;
    VB isBridge; VI visited, st, en, low;
    int clk = -1, edgeId = 0;
    Bridge(int emax, int vmax) : adj(vmax),
      isBridge(emax), visited(vmax), st(vmax),
      en(vmax), low(vmax) { ajs.reserve(emax); }
    void clearAll() {
      int n = adj.size();
      for(int i = 0; i < n; i++) {
        adj[i].clear(); visited[i] = st[i] = 0;
      }
      for(int i=0;i<edgeId;i++)isBridge[i]=false;
      clk = -1; edgeId = 0;
    }
    void dfs(int u, int parEdge) {
      visited[u] = 1; low[u] = st[u] = ++clk;
      for (int e : adj[u]) {
        if (e == parEdge) continue;
        int v = ajs[e].first ^ ajs[e].second ^ u;
        if (visited[v] == 1) {
          low[u] = min(low[u], st[v]);
        } else if(visited[v] == 0) {
          dfs(v, e); low[u] = min(low[u], low[v]);
        }
      }
      if(st[u] == low[u] && parEdge != -1) {
        isBridge[parEdge] = true;
      }
      en[u] = clk; visited[u] = 2;
    }
    void fndBriz() {
      int n = adj.size();
      for(int i = 0; i < n; i++) {
        if(visited[i] == 0) dfs(i, -1);
      }
    }
    bool isRplc(int eid, int u, int v) {
      if(!isBridge[eid]) return true;
      int a = ajs[eid].first, b = ajs[eid].second;
      if(st[a] > st[b]) swap(a, b);
      return (st[b] <= st[u] && st[u] <= en[b])
              != (st[b] <= st[v] && st[v] <= en[b]);
    }
    int addEdge(int u, int v) {
      ajs[edgeId] = {u, v};
      adj[u].emplace_back(edgeId);
      adj[v].emplace_back(edgeId);
      return edgeId++;
    }
  };
  vector< PII >ajs; int elmnts, grSz;
  GrafDual(int elmnts, int grSz) : bridge(elmnts,
    grSz), edgeMap(elmnts) {
    this->elmnts = elmnts; this->grSz = grSz;
    ajs.reserve(elmnts);
  }
  Bridge bridge; VI edgeMap;
  void updTknElmnt(const VB &tkn) {
    bridge.clearAll();
    for (int i = 0; i < elmnts; i++) if (!tkn[i]){
      edgeMap[i] = bridge.addEdge(ajs[i].first,
                                  ajs[i].second);
    }
    bridge.fndBriz();
  }
  bool canTkElmnt(int e) {
    return !bridge.isBridge[ edgeMap[e] ];
  }
  bool canEx(int rmv, int ins) {
    return bridge.isRplc(edgeMap[ins],
            ajs[rmv].first, ajs[rmv].second);
  }
};
bool augment(int elmnts, Mat *m1, Mat *m2,
    VB &tkn, const VB &source, const VB &sink) {
  VI parent(elmnts, -2), hidari, migi;
  hidari.reserve(elmnts); migi.reserve(elmnts);
  queue< int >q;
  for (int i = 0; i < elmnts; i++) {
    if (source[i]) {
      q.push(i); parent[i] = -1;
    }
    if (tkn[i]) hidari.push_back(i);
    else migi.push_back(i);
  }
  int connector = -1;
  while (!q.empty() && connector == -1) {
    int u = q.front(); q.pop();
    auto approach = [&](int v) {
      if (parent[v] == -2) {
        parent[v] = u;
        q.push(v);
        if (sink[v]) connector = v;
      }
    };
    if (tkn[u]) {
      for (int v : migi) if (m1->canEx(u, v)) {
        approach(v);
      }
    } else {
      for (int v : hidari) if (m2->canEx(v, u)) {
        approach(v);
      }
    }
  }
  if (connector == -1) return false;
  while (connector != -1) {
    tkn[connector] = tkn[connector] ^ 1;
    connector = parent[connector];
  }
  return true;
}
VB getBasisOfIntersection(int elmnts, Mat *m1,
    Mat *m2) {
  VB tkn(elmnts, false);
  while (true) {
    m1->updTknElmnt(tkn); m2->updTknElmnt(tkn);
    bool trivial=false,noSource=true,noSink=true;
    VB source(elmnts, false), sink(elmnts, false);
    for (int i = 0; i < elmnts; i++) {
      if (tkn[i]) continue;
      if (m1->canTkElmnt(i)) {
        source[i] = true; noSource = false;
      }
      if (m2->canTkElmnt(i)) {
        sink[i] = true; noSink = false;
      }
      if (source[i] && sink[i]) {
        tkn[i] = true; trivial = true; break;
      }
    }
    if(trivial)continue;if(noSource||noSink)break;
  if(!augment(elmnts,m1,m2,tkn,source,sink))break;
  }
  return tkn;
}
VI findEdgeDisjointSpanningTrees(const
  vector<PII>&ajs,int nodes,int trees) {
  int elmnts = ajs.size()*trees;
  GraphMat gm(elmnts, nodes*trees);
  ColorMat cm(elmnts, ajs.size());
  for (int i = 0; i < ajs.size(); i++) {
    PII p = ajs[i];
    for (int j = 0; j < trees; j++) {
     cm.elmntCol.push_back(i);gm.ajs.push_back(p);
     p.first += nodes; p.second += nodes;
    }
  }
  VB tkn = getBasisOfIntersection(elmnts,&gm,&cm);
  int on = 0; for (bool b : tkn) on += b;
  VI solution(ajs.size(), -1);
  if (on != trees*(nodes-1)) return solution;
  for (int i = 0; i < ajs.size(); i++) {
    for (int j = 0; j < trees; j++)
      if (tkn[i*trees+j]) {
      solution[i] = j;
    }
  }
  return solution;
}
bool weightedAugment(int elmnts, Mat *m1, Mat *m2,
  vector< CostType >costs, VB &tkn,
  const VB &source, const VB &sink) {
  VI parent(elmnts, -2), hidari, migi;
  hidari.reserve(elmnts); migi.reserve(elmnts);
  for (int i = 0; i < elmnts; i++) {
    if (tkn[i]) {
      hidari.push_back(i);
      costs[i] = -costs[i];
    } else {
      migi.push_back(i);
    }
  }
  vector< PII >exchangeajs;
  for (int u : hidari) {
    for (int v : migi) {
      if (m1->canEx(u, v)) {
        exchangeajs.emplace_back(u, v);
      }
      if (m2->canEx(u, v)) {
        exchangeajs.emplace_back(v, u);
      }
    }
  }
  vector< pair<CostType,int> >dist(elmnts,
    make_pair(INF, -1));
  for (int i = 0; i < elmnts; i++) if (source[i]){
    dist[i]=make_pair(costs[i], 0);parent[i] = -1;
  }
  for (int i = 0; i < elmnts; i++) {
    bool relaxed = false;
    for (PII p : exchangeajs) {
      if (parent[p.first] == -2) continue;
      pair< CostType, int >tmp = dist[p.first];
      tmp.first += costs[p.second]; tmp.second++;
      if (tmp < dist[p.second]) {
        relaxed = true; dist[p.second] = tmp;
        parent[p.second] = p.first;
      }
    }
    if (!relaxed) break;
  }
  int connector = -1;
  for (int i = 0; i < elmnts; i++)
    if (sink[i] && parent[i] != -2) {
    if (connector == -1||dist[i]<dist[connector]){
      connector = i;
    }
  }
  if (connector == -1) return false;
  while (connector != -1) {
    tkn[connector] = tkn[connector] ^ 1;
    connector = parent[connector];
  }
  return true;
}
/// returns rank+1 elmnts, minimum total costs for
/// a independent subset of size: 0, 1,.., rank
vector< CostType >weightedIntersection(int elmnts,
  Mat *m1, Mat *m2, const vector<CostType>&costs){
  VB tkn(elmnts, false);
  vector< CostType >minTotalCosts;
  while (true) {
    minTotalCosts.push_back(0);
    for (int i = 0; i < elmnts; i++) if (tkn[i]) {
      minTotalCosts.back() += costs[i];
    }
    m1->updTknElmnt(tkn); m2->updTknElmnt(tkn);
    bool noSource = true, noSink = true;
    VB source(elmnts, false), sink(elmnts, false);
    for (int i = 0; i < elmnts; i++) {
      if (tkn[i]) continue;
      if (m1->canTkElmnt(i)) {
        source[i] = true; noSource = false;
      }
      if (m2->canTkElmnt(i)) {
        sink[i] = true; noSink = false;
      }
    }
    if (noSource || noSink) break;
    if (!weightedAugment(elmnts, m1, m2, costs,
      tkn, source, sink)) break;
  }
  return minTotalCosts;
}
void solveURI_Honesty() {
  int n, m, k; int ti = 0;
  while (cin >> n >> m >> k) {
    GraphMat gm(m, n); ColorMat cm(m, k);
    for (int i = 0; i < m; i++) {
      int u, v, k; cin >> u >> v >> k;
      gm.ajs.emplace_back(u-1, v-1);
      cm.elmntCol.push_back(k-1);
    }
    VB basis = getBasisOfIntersection(m,&gm,&cm);
    int tkn = 0; for (bool b : basis) tkn += b;
    cout << "Instancia " << ++ti << "\n";
    if (tkn == n-1) cout << "sim" << "\n";
    else cout << "nao" << "\n";
    cout << "\n";
  }
}
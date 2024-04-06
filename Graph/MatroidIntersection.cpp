/** 
 * # Matroid Intersection - per increment O(r*n); r = rank of intersection, n = size of ground set
 * 1. Colorful Matroid : a set of elements is independent if each color's used frequency
 *  is not more than given respective threshold
 * 2. Graphic Matroid : a set of edges is independent if they form no cycle
 * 3. Binary Matroid : a set of elements is independent if it's linearly independent in GF(2)
 * 4. Graphic Dual Matroid : a set of edges is independent if complement edge set makes 
 *  the graph connected. # Caution: make sure empty set is independent for given input.
 * 
 * # Matroid Partitioning: findEdgeDisjointSpanningTrees() finds a given number of edge-disjoint 
 *   spanning trees in a graph. Similarly partitioning can be implemented for other types of matroids.
 * 
 * # Weighted Matroid Intersection: per increment O(r^2*n); in practice much faster because 
 *   of exchange graph being not very much customizable
 * 
 * # Notes:
 * => Everything 0-indexed
 * => During contest erase base Matroid to get better runtime (around 15% improvement)
 * => Tested on the problems specified here https://usaco.guide/adv/matroid-isect?lang=cpp
 * 
 * @author: solaimanope
*/
#include<bits/stdc++.h>
using namespace std;
typedef pair< int , int >PII;
typedef vector<int>VI;
typedef vector<bool>VB;
typedef long long CostType;
const CostType INF = 1e18;
const int BITSET_BITS = 60;
struct Graph {
  vector< vector< int > >edg;
  Graph(int nodes) : edg(nodes) { }
  void addEdge(int u, int v) {
    edg[u].push_back(v);}
  void clearGraph() {
    for (int i = 0; i < edg.size(); i++) 
     edg[i].clear();}
};
struct Matroid {
  virtual void updateTakenElements(const vector< bool > &taken) = 0;
  virtual bool canTakeElement(int e) = 0;
  virtual bool canExchange(int remove, int insert) = 0;
};
struct ColorfulMatroid : Matroid {
  vector< int >elementColor, canTakeAtMost, currentlyTaking;
  int elements, colors;
  ColorfulMatroid(int elements, int colors) : elements(elements), colors(colors),
      canTakeAtMost(colors, 1), currentlyTaking(colors, 0) {
    elementColor.reserve(elements);
  }
  void updateTakenElements(const vector< bool > &taken) {
    fill(currentlyTaking.begin(), currentlyTaking.end(), 0);
    for (int i = 0; i < elements; i++) if (taken[i]) {
      currentlyTaking[ elementColor[i] ]++;}
  }
  bool canTakeElement(int e) {
    int col = elementColor[e];
    return currentlyTaking[col] != canTakeAtMost[col];}
  bool canExchange(int remove, int insert){
    int colr = elementColor[remove];
    int coli = elementColor[insert];
    if (coli == colr) return true;
    return currentlyTaking[coli] != canTakeAtMost[coli];}
};
struct GraphicMatroid : Matroid {
  vector< PII >edges;
  int elements, graphSize;
  GraphicMatroid(int elements, int graphSize) : forest(graphSize) {
    this->elements = elements;
    this->graphSize = graphSize;
    edges.reserve(elements);
  }
  Graph forest;
  vector< int >start, finish, root;
  void treeDFS(int u, int p, int &tym) {
    start[u] = ++tym;
    for (int v : forest.edg[u]) {
      if (v == p) continue;
      root[v] = root[u];
      treeDFS(v, u, tym);
    }
    finish[u] = tym;
  }
  bool inSubtree(int u, int x) {
    return start[u] <= start[x] && finish[x] <= finish[u];
  }
  void updateTakenElements(const vector< bool > &taken) {
    forest.clearGraph();
    for (int i = 0; i < elements; i++) if (taken[i]) {
      forest.addEdge(edges[i].first, edges[i].second);
      forest.addEdge(edges[i].second, edges[i].first);
    }
    root = vector< int >(graphSize, -1);
    finish = start = vector< int >(graphSize, 0);
    int tym = -1;
    for (int i = 0; i < graphSize; i++) if (root[i] == -1) {
      root[i] = i;
      treeDFS(i, -1, tym);
    }
  }
  bool canTakeElement(int e) {
    return root[ edges[e].first ] != root[ edges[e].second ];
  }
  bool canExchange(int remove, int insert) {
    if (canTakeElement(insert)) return true;
    int u = edges[remove].first, p = edges[remove].second;
    if (start[p] > start[u]) swap(u, p);
    int x = edges[insert].first, y = edges[insert].second;
    return inSubtree(u, x) != inSubtree(u, y);
  }
};

struct BinaryMatroid : Matroid {
  typedef bitset< BITSET_BITS > BitSet;
  struct Basis {
    int usedBits;
    vector< BitSet >reduced;
    vector< BitSet >combination;
    Basis(int usedBits) : usedBits(usedBits), reduced(usedBits), combination(usedBits) {
      assert(BITSET_BITS == usedBits);
    }
    void clearAll() {
      for (int i = 0; i < usedBits; i++) {
        reduced[i].reset();
        combination[i].reset();
      }
    }
    BitSet canBeBuiltWith(BitSet x) {
      BitSet rt;
      for (int i = usedBits-1; i >= 0; i--) if (x.test(i)) {
        if (!reduced[i].test(i)) return BitSet();
        x ^= reduced[i];
        rt ^= combination[i];
      }
      return rt;
    }
    int addVector(BitSet x) {
      BitSet cm;
      for (int i = usedBits-1; i >= 0; i--) if (x.test(i)) {
        if (!reduced[i].test(i)) {
          reduced[i] = x;
          combination[i] = cm.set(i);
          return i;
        } else {
          x ^= reduced[i];
          cm ^= combination[i];
        }
      }
      return -1;
    }
  };
  vector< BitSet > rows;
  int elements, usedBits;
  BinaryMatroid(int elements, int usedBits) : elements(elements), usedBits(usedBits), 
      currentBasis(usedBits), cycle(elements), rowMap(elements) {
    rows.reserve(elements);
  }
  vector< BitSet >cycle;
  vector< int >rowMap;
  Basis currentBasis;
  void updateTakenElements(const vector< bool > &taken) {
    currentBasis.clearAll();
    for (int i = 0; i < elements; i++) if (taken[i]) {
      rowMap[i] = currentBasis.addVector(rows[i]);
    }
    for (int i = 0; i < elements; i++) if (!taken[i]) {
      cycle[i] = currentBasis.canBeBuiltWith(rows[i]);
    }
  }
  bool canTakeElement(int e) {
    return !cycle[e].any();
  }
  bool canExchange(int remove, int insert){
    if (canTakeElement(insert)) return true;
    return cycle[insert].test( rowMap[remove] );
  }
};

struct GraphicDual : Matroid {
  struct Bridge {
    vector< vector<int> >adj;
    vector< PII >edges;
    vector< bool >isBridge;
    vector< int >visited, st, en, low;
    int clk = -1, edgeId = 0;
    Bridge(int emax, int vmax) : adj(vmax), isBridge(emax), visited(vmax), 
        st(vmax), en(vmax), low(vmax) {
      edges.reserve(emax);
    }
    void clearAll() {
      int n = adj.size();
      for(int i = 0; i < n; i++) {
        adj[i].clear();
        visited[i] = st[i] = 0;
      }
      for(int i = 0; i < edgeId; i++) isBridge[i] = false;
      clk = -1;
      edgeId = 0;
    }
    void dfs(int u, int parEdge) {
      visited[u] = 1;
      low[u] = st[u] = ++clk;
      for (int e : adj[u]) {
        if (e == parEdge) continue;
        int v = edges[e].first ^ edges[e].second ^ u;
        if (visited[v] == 1) {
          low[u] = min(low[u], st[v]);
        } else if(visited[v] == 0) {
          dfs(v, e);
          low[u] = min(low[u], low[v]);
        }
      }
      if(st[u] == low[u] && parEdge != -1) {
        isBridge[parEdge] = true;
      }
      en[u] = clk;
      visited[u] = 2;
    }
    void findBridges() {
      int n = adj.size();
      for(int i = 0; i < n; i++) {
        if(visited[i] == 0) dfs(i, -1);
      }
    }
    bool isReplacable(int eid, int u, int v){
      if(!isBridge[eid]) return true;
      int a = edges[eid].first, b = edges[eid].second;
      if(st[a] > st[b]) swap(a, b);
      return (st[b] <= st[u] && st[u] <= en[b]) != (st[b] <= st[v] && st[v] <= en[b]);
    }
    int addEdge(int u, int v){
      edges[edgeId] = {u, v};
      adj[u].emplace_back(edgeId);
      adj[v].emplace_back(edgeId);
      return edgeId++;
    }
  };
  vector< PII >edges;
  int elements, graphSize;
  GraphicDual(int elements, int graphSize) : bridge(elements, graphSize), edgeMap(elements) {
    this->elements = elements;
    this->graphSize = graphSize;
    edges.reserve(elements);
  }
  Bridge bridge;
  vector< int >edgeMap;
  void updateTakenElements(const vector< bool > &taken) {
    bridge.clearAll();
    for (int i = 0; i < elements; i++) if (!taken[i]) {
      edgeMap[i] = bridge.addEdge(edges[i].first, edges[i].second);
    }
    bridge.findBridges();
  }
  bool canTakeElement(int e) {
    return !bridge.isBridge[ edgeMap[e] ];
  }
  bool canExchange(int remove, int insert) {
    return bridge.isReplacable(edgeMap[insert], edges[remove].first, edges[remove].second);
  }
};

bool augment(int elements, Matroid *m1, Matroid *m2, VB &taken, const VB &source, const VB &sink) {
  vector< int >parent(elements, -2), hidari, migi;
  hidari.reserve(elements);
  migi.reserve(elements);
  queue< int >q;
  for (int i = 0; i < elements; i++) {
    if (source[i]) {
      q.push(i);
      parent[i] = -1;
    }
    if (taken[i]) hidari.push_back(i);
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
    if (taken[u]) {
      for (int v : migi) if (m1->canExchange(u, v)) {
        approach(v);
      }
    } else {
      for (int v : hidari) if (m2->canExchange(v, u)) {
        approach(v);
      }
    }
  }
  if (connector == -1) return false;
  while (connector != -1) {
    taken[connector] = taken[connector] ^ 1;
    connector = parent[connector];
  }
  return true;
}

vector< bool >getBasisOfIntersection(int elements, Matroid *m1, Matroid *m2) {
  vector< bool >taken(elements, false);
  while (true) {
    m1->updateTakenElements(taken);
    m2->updateTakenElements(taken);
    bool trivial = false, noSource = true, noSink = true;
    vector< bool >source(elements, false), sink(elements, false);
    for (int i = 0; i < elements; i++) {
      if (taken[i]) continue;
      if (m1->canTakeElement(i)) {
        source[i] = true;
        noSource = false;
      }
      if (m2->canTakeElement(i)) {
        sink[i] = true;
        noSink = false;
      }
      if (source[i] && sink[i]) {
        taken[i] = true;
        trivial = true;
        break;
      }
    }
    if (trivial) continue;
    if (noSource || noSink) break;
    if (!augment(elements, m1, m2, taken, source, sink)) break;
  }
  return taken;
}

vector< int >findEdgeDisjointSpanningTrees(const vector< PII > &edges, int nodes, int trees) {
  int elements = edges.size()*trees;
  GraphicMatroid gm(elements, nodes*trees);
  ColorfulMatroid cm(elements, edges.size());
  for (int i = 0; i < edges.size(); i++) {
    PII p = edges[i];
    for (int j = 0; j < trees; j++) {
      cm.elementColor.push_back(i);
      gm.edges.push_back(p);
      p.first += nodes;
      p.second += nodes;
    }
  }
  vector< bool >taken = getBasisOfIntersection(elements, &gm, &cm);
  int on = 0;
  for (bool b : taken) on += b;

  vector< int >solution(edges.size(), -1);
  if (on != trees*(nodes-1)) return solution;
  for (int i = 0; i < edges.size(); i++) {
    for (int j = 0; j < trees; j++) if (taken[i*trees+j]) {
      solution[i] = j;
    }
  }
  return solution;
}

bool weightedAugment(int elements, Matroid *m1, Matroid *m2, vector< CostType >costs, VB &taken, const VB &source, const VB &sink) {
  vector< int >parent(elements, -2), hidari, migi;
  hidari.reserve(elements);
  migi.reserve(elements);
  for (int i = 0; i < elements; i++) {
    if (taken[i]) {
      hidari.push_back(i);
      costs[i] = -costs[i];
    } else {
      migi.push_back(i);
    }
  }
  vector< PII >exchangeEdges;
  for (int u : hidari) {
    for (int v : migi) {
      if (m1->canExchange(u, v)) {
        exchangeEdges.emplace_back(u, v);
      } 
      if (m2->canExchange(u, v)) {
        exchangeEdges.emplace_back(v, u);
      }
    }
  }
  vector< pair< CostType , int > >dist(elements, make_pair(INF, -1));
  for (int i = 0; i < elements; i++) if (source[i]) {
    dist[i] = make_pair(costs[i], 0);
    parent[i] = -1;
  }
  for (int i = 0; i < elements; i++) {
    bool relaxed = false;
    for (PII p : exchangeEdges) {
      if (parent[p.first] == -2) continue;
      pair< CostType , int >tmp = dist[p.first];
      tmp.first += costs[p.second];
      tmp.second++;
      if (tmp < dist[p.second]) {
        relaxed = true;
        dist[p.second] = tmp;
        parent[p.second] = p.first;
      }
    }
    if (!relaxed) break;
  }
  int connector = -1;
  for (int i = 0; i < elements; i++) if (sink[i] && parent[i] != -2) {
    if (connector == -1 || dist[i] < dist[connector]) {
      connector = i;
    }
  }
  if (connector == -1) return false;
  while (connector != -1) {
    taken[connector] = taken[connector] ^ 1;
    connector = parent[connector];
  }
  return true;
}

/// returns rank+1 elements, minimum total costs for a independent subset of size: 0, 1,.., rank
vector< CostType >weightedIntersection(int elements, Matroid *m1, Matroid *m2, const vector< CostType > &costs) {
  vector< bool >taken(elements, false);
  vector< CostType >minTotalCosts;
  while (true) {
    minTotalCosts.push_back(0);
    for (int i = 0; i < elements; i++) if (taken[i]) {
      minTotalCosts.back() += costs[i];
    }
    m1->updateTakenElements(taken);
    m2->updateTakenElements(taken);
    bool noSource = true, noSink = true;
    vector< bool >source(elements, false), sink(elements, false);
    for (int i = 0; i < elements; i++) {
      if (taken[i]) continue;
      if (m1->canTakeElement(i)) {
        source[i] = true;
        noSource = false;
      }
      if (m2->canTakeElement(i)) {
        sink[i] = true;
        noSink = false;
      }
    }
    if (noSource || noSink) break;
    if (!weightedAugment(elements, m1, m2, costs, taken, source, sink)) break;
  }
  return minTotalCosts;
}

void solveCodeChef_CNNCT2() {
  int t;
  cin >> t;

  while (t--) {
    int n, m;
    cin >> n >> m;

    GraphicMatroid gmA(m, n), gmB(m, n);
    for (int i = 0; i < m; i++) {
      int u, v;
      cin >> u >> v;
      gmA.edges.emplace_back(u-1, v-1);
    }
    for (int i = 0; i < m; i++) {
      int u, v;
      cin >> u >> v;
      gmB.edges.emplace_back(u-1, v-1);
    }

    vector< bool >base = getBasisOfIntersection(m, &gmA, &gmB);
    int ans = 2*(n-1);
    for (bool x : base) ans -= x;
    cout << ans << "\n";
  }
}

void solveCodeChef_CNNCT2__dual() {
  int t;
  cin >> t;

  while (t--) {
    int n, m;
    cin >> n >> m;

    GraphicDual dgmA(m, n), dgmB(m, n);
    for (int i = 0; i < m; i++) {
      int u, v;
      cin >> u >> v;
      dgmA.edges.emplace_back(u-1, v-1);
    }
    for (int i = 0; i < m; i++) {
      int u, v;
      cin >> u >> v;
      dgmB.edges.emplace_back(u-1, v-1);
    }

    vector< bool >base = getBasisOfIntersection(m, &dgmA, &dgmB);
    int ans = m;
    for (bool x : base) ans -= x;
    cout << ans << "\n";
  }
}

void solveURI_Honesty() {
  int n, m, k;
  int ti = 0;
  while (cin >> n >> m >> k) {
    GraphicMatroid gm(m, n);
    ColorfulMatroid cm(m, k);
    for (int i = 0; i < m; i++) {
      int u, v, k;
      cin >> u >> v >> k;
      gm.edges.emplace_back(u-1, v-1);
      cm.elementColor.push_back(k-1);
    }

    vector< bool >basis = getBasisOfIntersection(m, &gm, &cm);
    int taken = 0;
    for (bool b : basis) taken += b;

    cout << "Instancia " << ++ti << "\n";
    if (taken == n-1) cout << "sim" << "\n";
    else cout << "nao" << "\n";
    cout << "\n";
  }
}

void solveSPOJ_COIN() {
  while (true) {
    int r;
    cin >> r;
    if (r == 0) break;
    vector< int >w;
    vector< array< int , 4 > >rounds;    
    for (int i = 0; i < r; i++) {
      rounds.emplace_back();
      for (int &x : rounds.back()) {
        cin >> x;
        w.push_back(x);
      }
    }
    sort(w.begin(), w.end());
    w.erase(unique(w.begin(), w.end()), w.end());

    GraphicMatroid gm(r+r, w.size());
    ColorfulMatroid cm(r+r, r);
    for (int i = 0; i < r; i++) {
      for (int &x : rounds[i]) {
        x = lower_bound(w.begin(), w.end(), x) - w.begin();
      }
      gm.edges.emplace_back(rounds[i][0], rounds[i][1]);
      gm.edges.emplace_back(rounds[i][2], rounds[i][3]);
      cm.elementColor.push_back(i);
      cm.elementColor.push_back(i);
    }

    vector< bool >basis = getBasisOfIntersection(r+r, &gm, &cm);

    int ans = 0;
    for (bool b : basis) ans += b;
    cout << ans*2 << "\n";
  }
}

void solveCodeForces_Seollal() {
  struct Info {
    bool hori;
    int u, v;
    int black;
  };
  int dx[] = {1, 0, -1, 0};
  int dy[] = {0, 1, 0, -1};

  int t;
  cin >> t;

  while (t--) {
    int n, m;
    cin >> n >> m;
    auto encode = [&](int i, int j) {
      return i*m+j;
    };
    auto decode = [&](int x) {
      return PII(x/m, x%m);
    };

    vector< int >fb;
    vector< string >grid(n);
    int freeCells = 0;
    for (int i = 0; i < n; i++) {
      cin >> grid[i];
      for (int j = i&1; j < m; j += 2) {
        if (grid[i][j] == 'X') continue;
        fb.push_back(encode(i, j));
      }
      for (char c : grid[i]) freeCells += c == 'O';
    }

    vector< int >options;
    for (int x : fb) {
      PII p = decode(x);
      options.push_back(0);
      for (int k = 0; k < 4; k++) {
        int i = p.first+dx[k];
        int j = p.second+dy[k];
        if (i < 0 || i >= n || j < 0 || j >= m) continue;
        if (grid[i][j] == 'O') options.back()++;
      }
    }
    int idx = 0;
    for (int i = 1; i < options.size(); i++) {
      if (options[i] <= options[idx]) idx = i;
    }
    if (idx > 0 && options[idx] < 2) {
      cout << "NO" << "\n";
      continue;
    }

    vector< Info >edges;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j+1 < m; j++) {
        if (grid[i][j] == 'O' && grid[i][j+1] == 'O') {
          int black = lower_bound(fb.begin(), fb.end(), encode(i, j+((i^j)&1))) - fb.begin();
          edges.push_back({true, encode(i, j), encode(i, j+1), black});
        }
      }
    }
    for (int i = 0; i+1 < n; i++) {
      for (int j = 0; j < m; j++) {
        if (grid[i][j] == 'O' && grid[i+1][j] == 'O') {
          int black = lower_bound(fb.begin(), fb.end(), encode(i+((i^j)&1), j)) - fb.begin();
          edges.push_back({false, encode(i, j), encode(i+1, j), black});
        }
      }
    }
    GraphicDual dgm(edges.size(), n*m);
    ColorfulMatroid cm(edges.size(), fb.size());
    for (int i = 0; i < edges.size(); i++) {
      dgm.edges.emplace_back(edges[i].u, edges[i].v);
      cm.elementColor.push_back(edges[i].black);
    }
    for (int i = 1; i < fb.size(); i++) {
      cm.canTakeAtMost[i] = options[i] - 2;
    }
    cm.canTakeAtMost[0] = 2;

    vector< bool >basis = getBasisOfIntersection(edges.size(), &dgm, &cm);
    int off = 0;
    for (bool b : basis) off += !b;

    if (off != freeCells - 1) {
      cout << "NO" << "\n";
      continue;
    }

    vector< string >soln(n+n-1, string(m+m-1, ' '));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        soln[i*2][j*2] = grid[i][j];
      }
    }
    for (int i = 0; i < edges.size(); i++) {
      if (basis[i]) continue;
      PII p = decode(edges[i].u);
      if (edges[i].hori) {
        soln[p.first*2][p.second*2+1] = '-';
      } else {
        soln[p.first*2+1][p.second*2] = '|';
      }
    }

    cout << "YES\n";
    for (string &s : soln) cout << s << "\n";
  }
}

void solveGym_PickYourOwnNim() {
  typedef long long LL;
  vector< pair< int , LL > >options;
  int n;
  cin >> n;
  for (int i = 0; i < n; i++) {
    LL x;
    cin >> x;
    options.emplace_back(i, x);
  }
  int m;
  cin >> m;
  for (int i = 0; i < m; i++) {
    int k;
    cin >> k;
    BinaryMatroid::Basis basis(BITSET_BITS);
    while (k--) {
      LL x;
      cin >> x;
      if (basis.addVector(x) != -1) {
        options.emplace_back(n+i, x);
      }
    }
  }

  BinaryMatroid bm(options.size(), BITSET_BITS);
  ColorfulMatroid cm(options.size(), n+m);
  for (int i = 0; i < options.size(); i++) {
    cm.elementColor.push_back(options[i].first);
    bm.rows.push_back(options[i].second);
  }

  vector< bool >taken = getBasisOfIntersection(options.size(), &bm, &cm);

  int on = 0;
  for (bool b : taken) on += b;
  if (on != n+m) {
    cout << -1 << "\n";
    return;
  }
  for (int i = n; i < options.size(); i++) {
    if (taken[i]) {
      cout << options[i].second << "\n";
    }
  }
}

void solveUVa_Ambiguous_Forests() {
  typedef array< int , 4 > Quad;
  int t;
  cin >> t;

  while (t--) {
    int n1, n2, m;
    cin >> n1 >> n2 >> m;

    auto augment = [](int elements, Matroid *m1, Matroid *m2, VB &taken, const VB &source, const VB &sink, int exclude) {
      vector< int >parent(elements, -2), hidari, migi;
      hidari.reserve(elements);
      migi.reserve(elements);
      queue< int >q;
      for (int i = 0; i < elements; i++) {
        if (source[i]) {
          q.push(i);
          parent[i] = -1;
        }
        if (!taken[i]) migi.push_back(i);
        else if (i < exclude) hidari.push_back(i);
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
        if (taken[u]) {
          for (int v : migi) if (m1->canExchange(u, v)) {
            approach(v);
          }
        } else {
          for (int v : hidari) if (m2->canExchange(v, u)) {
            approach(v);
          }
        }
      }
      if (connector == -1) return false;
      while (connector != -1) {
        taken[connector] = taken[connector] ^ 1;
        connector = parent[connector];
      }
      return true;
    };

    auto canAddMore = [&](vector< Quad >options, vector< Quad > must) {
      int elements = options.size() + must.size();
      GraphicMatroid gm1(elements, n1);
      GraphicMatroid gm2(elements, n2);

      vector< bool >taken;
      for (Quad q : options) {
        gm1.edges.emplace_back(q[0], q[1]);
        gm2.edges.emplace_back(q[2], q[3]);
        taken.push_back(false);
      }
      for (Quad q : must) {
        gm1.edges.emplace_back(q[0], q[1]);
        gm2.edges.emplace_back(q[2], q[3]);
        taken.push_back(true);
      }
      while (true) {
        gm1.updateTakenElements(taken);
        gm2.updateTakenElements(taken);
        bool trivial = false, noSource = true, noSink = true;
        vector< bool >source(elements, false), sink(elements, false);
        for (int i = 0; i < elements; i++) {
          if (taken[i]) continue;
          if (gm1.canTakeElement(i)) {
            source[i] = true;
            noSource = false;
          }
          if (gm2.canTakeElement(i)) {
            sink[i] = true;
            noSink = false;
          }
          if (source[i] && sink[i]) {
            taken[i] = true;
            trivial = true;
            break;
          }
        }
        if (trivial) continue;
        if (noSource || noSink) break;
        if (!augment(elements, &gm1, &gm2, taken, source, sink, options.size())) break;
      }
      int rt = 0;
      for (int i = 0; i < options.size(); i++) rt += taken[i];
      return rt;
    };

    vector< Quad >options(m);
    for (int i = 0; i < m; i++) {
      for (int &x : options[i]) cin >> x;
    }
    vector< Quad > must;
    int ans = canAddMore(options, must);
    // cout << "ans = " << ans << endl;
    vector< int >soln(ans, -1);
    for (int i = m-1, run = ans; i >= 0 && run > 0; i--) {
      Quad tmp = options.back();
      options.pop_back();
      int without = canAddMore(options, must);
      // cout << "without " << i << " -> " << without << endl;
      if (without < run) {
        must.push_back(tmp);
        soln[--run] = i;
        // cout << "soln[" << run << "] = " << i << endl;
      }
    }
    cout << ans;
    for (int x : soln) cout << " " << x;
    cout << "\n";
  }
}

void solveCodeJam_DatacenterDuplex() {
  int dx[] = {1, 0, -1, 0};
  int dy[] = {0, 1, 0, -1};

  int t;
  cin >> t;

  for (int ti = 1; ti <= t; ti++) {
    int R, C;
    cin >> R >> C;
    vector< string >grid(R);
    for (string &s : grid) cin >> s;

    vector< vector< int > >visited(R, vector< int >(C, -1));

    function< void(int, int) >dfs = [&](int i, int j) {
      for (int k = 0; k < 4; k++) {
        int x = i+dx[k];
        int y = j+dy[k];
        if (x < 0 || y < 0 || x >= R || y >= C) continue;
        if (grid[x][y] == grid[i][j] && visited[x][y] == -1) {
          visited[x][y] = visited[i][j];
          dfs(x, y);
        }
      }
    };

    int comps = 0;
    for (int i = 0; i < R; i++) {
      for (int j = 0; j < C; j++) {
        if (visited[i][j] == -1) {
          visited[i][j] = comps++;
          dfs(i, j);
        }
      }
    }

    struct Info {
      bool anti;
      int i, j;
      int u, v;
      int c;
    };

    vector< Info >options;
    int colors = 0;
    for (int i = 0; i+1 < R; i++) {
      for (int j = 0; j+1 < C; j++) {
        bool pushed = false;
        if (grid[i][j] == grid[i+1][j+1] && visited[i][j] != visited[i+1][j+1]) {
          pushed = true;
          options.push_back({false, i, j, visited[i][j], visited[i+1][j+1], colors});
        }
        if (grid[i][j+1] == grid[i+1][j] && visited[i][j+1] != visited[i+1][j]) {
          pushed = true;
          options.push_back({true, i, j, visited[i][j+1], visited[i+1][j], colors});
        }
        if (pushed) colors++;
      }
    }

    ColorfulMatroid cm(options.size(), colors);
    GraphicMatroid gm(options.size(), comps);
    for (Info i : options) {
      gm.edges.emplace_back(i.u, i.v);
      cm.elementColor.push_back(i.c);
    }
    
    vector< bool >taken = getBasisOfIntersection(options.size(), &gm, &cm);
    int on = 0;
    for (bool b : taken) on += b;
    
    if (on < comps-2) {
      cout << "Case #" << ti << ": IMPOSSIBLE" << "\n";
    } else {
      cout << "Case #" << ti << ": POSSIBLE" << "\n";
      vector< string >soln(R-1, string(C-1, '.'));
      for (int i = 0; i < options.size(); i++) if (taken[i]) {
        if (options[i].anti) {
          soln[ options[i].i ][ options[i].j ] = '/';
        } else {
          soln[ options[i].i ][ options[i].j ] = '\\';
        }
      }
      for (string s : soln) cout << s << "\n";
    }
  }
}

void solve_dmopc19c3p6() {
  int n, m;
  cin >> n >> m;
  vector< PII >edges(m);
  for (PII &p : edges) {
    cin >> p.first >> p.second;
    p.first--;
    p.second--;
  }
  vector< int >solution = findEdgeDisjointSpanningTrees(edges, n, 3);
  if (*max_element(solution.begin(), solution.end()) == -1) {
    cout << -1 << "\n";
    return;
  }
  for (int x : solution) cout << x+1;
  cout << "\n";
}

void solveKattis_rainbowgraph() {
  int n, m;
  cin >> n >> m;

  GraphicDual red(m, n), blue(m, n);
  vector< CostType >costs;
  CostType total = 0;
  for (int i = 0; i < m; i++) {
    int a, b, w; char c;
    cin >> a >> b >> w >> c;
    a--; b--;
    if (c == 'R') {
      red.edges.emplace_back(a, b);
      blue.edges.emplace_back(0, 0);
    } else if (c == 'B') {
      blue.edges.emplace_back(a, b);
      red.edges.emplace_back(0, 0);
    } else {
      red.edges.emplace_back(a, b);
      blue.edges.emplace_back(a, b);
    }
    costs.push_back(-w);
    total += w;
  }

  Graph gr(n);
  vector< bool >visited(n, false);
  function<int(int)>compSize = [&](int u) {
    int ans = 1;
    visited[u] = true;
    for (int v : gr.edg[u]) if (!visited[v]) {
      ans += compSize(v);
    }
    return ans;
  };

  for (PII p : red.edges) {
    gr.addEdge(p.first, p.second);
    gr.addEdge(p.second, p.first);
  }
  if (compSize(0) < n) {
    for (int i = 0; i < m; i++) cout << -1 << "\n";
    return;
  }
  gr.clearGraph();
  fill(visited.begin(), visited.end(), false);
  for (PII p : blue.edges) {
    gr.addEdge(p.first, p.second);
    gr.addEdge(p.second, p.first);
  }
  if (compSize(0) < n) {
    for (int i = 0; i < m; i++) cout << -1 << "\n";
    return;
  }

  vector< CostType >soln = weightedIntersection(m, &red, &blue, costs);
  
  for (int i = m; i > 0; i--) {
    if (soln.size() < i) cout << -1 << "\n";
    else cout << total+soln[i-1] << "\n";
  }
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);

  //solveCodeChef_CNNCT2();
  //solveCodeChef_CNNCT2__dual();
  //solveURI_Honesty();
  //solveSPOJ_COIN();
  //solveCodeForces_Seollal();
  //solveGym_PickYourOwnNim();
  //solveUVa_Ambiguous_Forests();
  //solveCodeJam_DatacenterDuplex();
  //solve_dmopc19c3p6();
  solveKattis_rainbowgraph();

  return 0;
}

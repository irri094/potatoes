const int vmax = 2e5+10, emax = 2e5+10;
namespace Bridge {///edge, nodes, comps 1 indexed
  vector<int> adj[vmax]; /// edge-id
  pair<int, int> edges[emax];   /// (u, v)
  bool isBridge[emax];
  int visited[vmax];  ///0-unvis,1-vising,2-vis
  int st[vmax], low[vmax], clk = 0, edgeId = 0;
  /// For bridge tree components
  int who[vmax], compId = 0;
  vector<int> stk;
  /// For extra end time calc
  int en[vmax];
  void dfs(int u, int parEdge) {
    visited[u] = 1; low[u] = st[u] = ++clk;
    stk.push_back(u);
    for (auto e : adj[u]) {
      if (e == parEdge) continue;
      int v=edges[e].first^edges[e].second^u;
      if (visited[v] == 1) {
        low[u] = min(low[u], st[v]);
      } else if(visited[v] == 0){
        dfs(v, e); low[u] = min(low[u], low[v]);
      }
    }
    visited[u] = 2;
    if(st[u] == low[u]){/// found
      ++compId; int cur;
      do{
        cur = stk.back(); stk.pop_back();
        who[cur] = compId;
      }while(cur != u);
      if(parEdge != -1){isBridge[parEdge] = true;}
    }
    en[u] = clk;
  }
  void clearAll(int n){
    for(int i = 0; i<=n; i++) {
      adj[i].clear(); visited[i] = st[i] = 0; }
    for(int i = 0; i<=edgeId; i++) isBridge[i]=0;
    clk = compId = edgeId = 0;
  }
  void findBridges(int n){
    for(int i = 1; i<=n; i++){
      if(visited[i] == 0) dfs(i, -1); }
  }
  bool isReplacable(int eid, int u, int v){
    if(!isBridge[eid]) return true;
    int a=edges[eid].first,b=edges[eid].second;
    if(st[a] > st[b]) swap(a, b);
    return (st[b] <= st[u] && st[u] <= en[b])
    != (st[b] <= st[v] && st[v] <= en[b]);
  }
  void addEdge(int u, int v){
    edgeId++; edges[edgeId] = {u, v};
    adj[u].emplace_back(edgeId);
    adj[v].emplace_back(edgeId);
  }
}

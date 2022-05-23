const int vmax = 1e4+10, emax = 1e5+10;
namespace Cutpoint { /// For BCTree, no self edge
/// edge, nodes, components 1-indexed
  vector<int> adj[vmax];   ///edge-id
  pair<int, int> edges[emax];   /// (u, v)
  bool isCutpoint[vmax];
  int visited[vmax];///0-unvis, 1-vising, 2-vis
  int st[vmax], low[vmax], clk = 0, edgeId = 0;
  /// For block components (i.e. edges)
  int who[emax], compId = 0;
  vector<int> stk;
  /// For extra end time calc
  int en[vmax];
  void dfs(int u, int parEdge) {
    visited[u] = 1; low[u] = st[u] = ++clk;
    int ch_cnt = (parEdge != -1);
    for (auto e : adj[u]) {
      if (e == parEdge) continue;
      int v = edges[e].first^edges[e].second^u;
      if (visited[v] == 1) {
        stk.push_back(e);
        low[u] = min(low[u], st[v]);
      } else if(visited[v] == 0){
        stk.push_back(e); dfs(v, e);
        low[u] = min(low[u], low[v]);
        if(low[v] >= st[u]){
          ++ch_cnt; ++compId;
          int cur;
          do{
            cur = stk.back(); stk.pop_back();
            who[cur] = compId;
          }while(cur != e);
        }
      }
    } visited[u] = 2;
    if(ch_cnt > 1){ isCutpoint[u] = true;}
    en[u] = clk;
  }
  void clearAll(int n){
    for(int i = 0; i<=n; i++) {
      adj[i].clear(); visited[i] = st[i] = 0; }
    for(int i=0;i<=n;i++)isCutpoint[i]=false;
    clk = compId = edgeId = 0;
  }
  void findCutpoints(int n){
    for(int i = 1; i<=n; i++){
      if(visited[i] == 0) dfs(i, -1); }
  }
  void addEdge(int u, int v){
    edgeId++; edges[edgeId] = {u, v};
    adj[u].emplace_back(edgeId);
    adj[v].emplace_back(edgeId);
  }
}

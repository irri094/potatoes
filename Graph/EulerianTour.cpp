const int vmax = 1e5+10, emax = 2e5+10;
namespace Euler { ///nodes, edges 1 indexed [1, n]
  /// call clear(vertex count) to clear stuff
  pair<int, int> edges[emax];
  bool used[emax]; /// used edges
  int ecnt = 0;
  vector<int> adj[vmax];
  int ptr[vmax];  /// curr pointer at adjlist
  void addEdge(int u, int v, const bool&
    directed = false){
    edges[++ecnt] = {u, v};
    adj[u].push_back(ecnt);
    if (!directed)  adj[v].push_back(ecnt);
  }
  vector<pair<int, int>> stk; ///(node, edge)
  vector<int>eulerTour(int st=edges[ecnt].first){
    vector<int> tour;
    if(ecnt == 0) return tour;
    stk.emplace_back(st, -1);
    while(!stk.empty()){
      int u = stk.back().first;
      for(int &i = ptr[u]; i<adj[u].size(); i++){
        int e = adj[u][i];
        if(used[e]) continue;
        used[e] = true;
        int v = u^edges[e].first^edges[e].second;
        stk.emplace_back(v, e); break;
      }
      if(ptr[u] == adj[u].size()){
        tour.push_back(stk.back().second);
        stk.pop_back();
      }
    }
    tour.pop_back();
    reverse(tour.begin(), tour.end());
    return tour;
  }
  vector<int> eulerPath(int st, int en){
    addEdge(en, st, true);
    auto path = eulerTour(st);
    if(ecnt != path.back()) {
     auto it=find(path.begin(),path.end(),ecnt)+1;
      rotate(path.begin(), it, path.end());
    }
    path.pop_back(); --ecnt;
    adj[en].pop_back(); adj[st].pop_back();
    return path;
  }
  void clear(int n){
    for(int i=0;i<=n;i++)adj[i].clear(),ptr[i]=0;
    for(int i = 0; i<=ecnt; i++) used[i] = false;
    ecnt = 0;
  }
}

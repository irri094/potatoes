namespace Dinic {     typedef long long LL;
const int N = 5005, K = 60; const LL INF = 1e18;
struct Edge { int frm, to; LL cap, flow; };
int s, t, n, level[N], ptr[N];
vector<Edge> edges; vector<int> adj[N];
void init(int nodes) {
  n = nodes; edges.clear();
  for (int i=0; i<n; i++)  adj[i].clear();
}
///adding undirected edge call addEdge(u,v,c,c);
int addEdge(int a, int b, LL cap, LL revcap = 0) {
  edges.push_back({a, b, cap, 0});
  edges.push_back({b, a, revcap, 0});
  adj[a].push_back(edges.size()-2);
  adj[b].push_back(edges.size()-1);
  return edges.size()-2;
}
bool bfs(LL lim) {
  fill(level, level+n, -1); level[s] = 0;
  queue<int> q; q.push(s);
  while (!q.empty() && level[t] == -1) {
    int v = q.front(); q.pop();
    for (int id: adj[v]) {
      Edge e = edges[id];
      if (level[e.to]==-1 && e.cap-e.flow>=lim) {
        q.push(e.to); level[e.to] = level[v] + 1;
      }
    }
  }
  return level[t] != -1;
}
LL dfs(int v, LL flow) {
  if (v == t || !flow)    return flow;
  for (; ptr[v] < adj[v].size(); ptr[v]++) {
    int eid = adj[v][ptr[v]];
    Edge &e = edges[eid];
    if (level[e.to] != level[v] + 1)  continue;
   if(LL pushed=dfs(e.to,min(flow,e.cap-e.flow))){
      e.flow+=pushed; edges[eid^1].flow -= pushed;
      return pushed;
    }
  }
  return 0;
}
LL maxFlow(int source,int sink,bool SCALING=false){
  s = source, t = sink;
  long long flow = 0;
  for (LL lim=SCALING?(1LL<<K):1; lim>0; lim>>=1){
    while (bfs(lim)) {
      fill(ptr, ptr+n, 0);
      while (LL pushed = dfs(s,INF)) flow+=pushed;
    }
  }
  return flow;
}
bool leftOfMinCut(int x) {return level[x] != -1;}
vector<vector<LL>> allPairFlow(vector<Edge>&tree){
  tree.clear(); vector<int> par(n);
  vector<vector<LL>> flow(n, vector<LL>(n, INF));
  for (int i=1; i<n; i++) {
    for (auto &e: edges)  e.flow = 0;
    LL f = maxFlow(i, par[i]);
    tree.push_back({i, par[i], f});
    for (int j=i+1; j<n; j++)
      if(par[j]==par[i]&&leftOfMinCut(j))par[j]=i;
    flow[i][par[i]] = flow[par[i]][i] = f;
    for (int j=0; j<i; j++)
      if (j != par[i])  flow[i][j]=flow[j][i]
                         =min(f,flow[par[i]][j]);
  }
  return flow;
} }

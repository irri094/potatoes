namespace SCC { //Everything 0-indexed.
const int N = 2e6+7; int which[N], vis[N], cc;
vector<int> adj[N], adjr[N]; vector<int> order;
void addEdge(int u, int v) {
  adj[u].push_back(v); adjr[v].push_back(u);
}
void dfs1(int u){
  if (vis[u]) return; vis[u] = true;
  for(int v: adj[u]) dfs1(v); order.push_back(u);
}
void dfs2(int u, int id) {
  if(vis[u]) return; vis[u] = true;
  for(int v: adjr[u])  dfs2(v, id); which[u] = id;
}
int last = 0;
void findSCC(int n) {
  cc=0,last=n; order.clear(); fill(vis, vis+n, 0);
  for(int i=0; i<n; i++) if(!vis[i])   dfs1(i);
  reverse(order.begin(), order.end());
  fill(vis, vis+n, 0);
  for (int u: order) {
    if (vis[u])   continue; dfs2(u, cc); ++cc;
  }
}
void clear() {
  for (int i=0; i<last; i++)
    adj[i].clear(), adjr[i].clear();
} }
struct TwoSat {
  int n; int vars = 0; vector<bool> ans;
  TwoSat(int n) : n(n), ans(n) {
    SCC::clear(); vars = 2*n;
  }
  void implies(int x, int y) {
    SCC::addEdge(x, y); SCC::addEdge(y^1, x^1);
  }
  void OR(int x, int y) {
    SCC::addEdge(x^1, y); SCC::addEdge(y^1, x);
  }
  void XOR(int x, int y) {
    implies(x, y^1); implies(x^1, y);
  }
  void atmostOne(vector<int> v) {
    int k = v.size();
    for (int i=0; i<k; i++) {
      if (i+1<k)    implies(vars+2*i, vars+2*i+2);
      implies(v[i], vars+2*i);
      if (i>0)    implies(v[i], vars+2*i-1);
    }
    vars += 2*k;
  }
  bool solve() {
    SCC::findSCC(vars); ans.resize(vars/2);
    for (int i=0; i<vars; i+=2) {
      if (SCC::which[i]==SCC::which[i+1])return 0;
      if (i<2*n)
        ans[i/2] = SCC::which[i]>SCC::which[i+1];
    }
    return true;
  }
};

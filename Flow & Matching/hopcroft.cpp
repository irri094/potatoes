namespace HopcroftKarp { // 1-indexed. L,R indep
const int maxN=1e5+7,maxM=1e5+7; int n, m, match;
int vis[maxN], level[maxN], ml[maxN], mr[maxM];
vector<int> edge[maxN];
void init(int N, int M) { //N=left nodes, M = right
  n = N, m = M;
  for (int i=1;i<=n;i++) edge[i].clear(),ml[i]=-1;
  for (int i=1;i<=m;i++)  mr[i] = -1;
  match = 0;
}
void add(int u, int v) {
  edge[u].push_back(v);
  if(ml[u]==-1&&mr[v]==-1)ml[u]=v,mr[v]=u,match++;
}
bool dfs(int u) {
  vis[u] = true;
  for (int x: edge[u]) {
    int v = mr[x];
    if (v == -1 || (!vis[v]
              && level[u] < level[v] && dfs(v))) {
      ml[u] = x; mr[x] = u; return true;
    }
  }
  return false;
}
int matching() {
  while (true) {
    queue<int> q;
    for (int i = 1; i <= n; ++i) {
      if (ml[i] == -1)  level[i] = 0, q.push(i);
      else        level[i] = -1;
    }
    while (!q.empty()) {
      int u = q.front(); q.pop();
      for (int x: edge[u]) {
        int v = mr[x];
        if (v != -1 && level[v] < 0) {
          level[v] = level[u] + 1; q.push(v);
        }
      }
    }
    for (int i = 1; i <= n; ++i) vis[i] = false;
    int d = 0;
    for(int i=1;i<=n;++i)if(ml[i]==-1&&dfs(i))++d;
    if (d == 0) return match; match += d;
  }
} }

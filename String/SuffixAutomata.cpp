namespace Automata {
const int N = 1e6+7, K = 26;
int len[2*N], link[2*N], nxt[2*N][K], sz, last;
void init(int n) {
  for (int i=0; i<=2*n; i++)
    fill(nxt[i], nxt[i]+K, -1);
  len[0] = 0; link[0] = -1; sz = 1; last = 0;
}
void add(char ch) {
  int c = ch-'a', cur = sz++; //create new node
  len[cur] = len[last]+1;
  int u = last;
  while (u != -1 && nxt[u][c] == -1) {
    nxt[u][c] = cur; u = link[u]; }

  if (u == -1) { link[cur] = 0; }
  else {
    int v = nxt[u][c];
    if (len[v] == len[u]+1) { link[cur] = v; }
    else {
      int clone = sz++;  //create node by cloning
      len[clone] = 1+len[u]; link[clone]=link[v];
      for (int i=0; i<K; i++)
        nxt[clone][i] = nxt[v][i];
      while (u != -1 && nxt[u][c] == v) {
        nxt[u][c] = clone; u = link[u];
      }
      link[v] = link[cur] = clone;
    }
  } last = cur;
}
vector<int> edge[2*N];
void makeEdge() { ///Make Suffix Link Tree edges
  for (int i=0; i<sz; i++) {
    edge[i].clear();
    for (int j=0; j<K; j++)
      if (nxt[i][j]!=-1)    edge[i].push_back(j);
  }
} }

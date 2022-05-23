/// O(n^3) OUTPUT: (min cut value, nodes in left)
pair<LL,vector<int>>glMinCut(vector<vector<LL>>c){
  int N = c.size();
  LL best_w = -1;
  vector<int> used(N), cut, best_cut;
  for (int phase = N-1; phase >= 0; phase--) {
    vector<LL> w = c[0];
    vector<int> vis = used;
    int prev, last = 0;
    for (int i = 0; i < phase; i++) {
      prev = last;
      last = -1;
      for (int j = 1; j < N; j++)
        if (!vis[j] && (last==-1 || w[j]>w[last]))
          last = j;
      if (i == phase-1) {
       for(int j=0;j<N;j++)c[prev][j]+=c[last][j];
       for(int j=0;j<N;j++)c[j][prev] =c[prev][j];
        used[last] = true; cut.push_back(last);
        if (best_w==-1 || w[last] < best_w) {
          best_cut = cut; best_w = w[last];
        }
      } else {
        for (int j=0; j<N; j++) w[j]+=c[last][j];
        vis[last] = true;
      }
    }
  } return make_pair(best_w, best_cut);
}

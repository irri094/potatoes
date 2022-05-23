const double EPS = 1e-9;
typedef vector<vector<double>> VVD;
int gauss(VVD ar, int e, VVD &res){
  int n = ar.size(), m = ar[0].size()-1;
  vector<int> pos(m, -1);
  for (int j=0, i=0; j<m && i<n; j++) {
    int p = i;
    for (int k=i; k<n; k++){
      if (abs(ar[k][j]) > abs(ar[p][j])) p = k;
    }
    if (abs(ar[p][j]) > EPS){
      pos[j] = i; swap(ar[p], ar[i]);
      for (int k=0; k<n; k++){
        if (k == i) continue;
        double x = ar[k][j]/ar[i][j];
        for(int l=j;l<m+e;l++)ar[k][l]-=ar[i][l]*x;
      } i++;
    }
  }
  int free_var = 0;
  for (int i=0;i<m; i++) free_var += (pos[i]==-1);
  for (int k=0; k<e; k++) {
    res.emplace_back(vector<double>(m));
    for (int i=0; i<m; i++)
      res.back()[i]=ar[pos[i]][m+k]/ar[pos[i]][i];
//    for (int i = 0; i < n; i++) {
//      double val = 0;
//      for (int j=0; j<m; j++)
//        val += res.back()[j]*ar[i][j];
//   if (abs(val-ar[i][m])>EPS)res.back().clear();
//  }
  } return free_var;
}

//vector<pair<int,ll> info;
//(original row, mask of ranked current rows)
int xorGaussian(vector<ll> &rows){
  int r = 0, n = (int)rows.size();
  //info.resize(n);
  //for(int i = 0; i<n; i++) info[i] = {i, 0ll};
  for(int c = 63; c>=0; c--){
    int p;
    for(p = r; p<n; p++) if(rows[p]&bt(c)) break;
    if(p == n) continue;
    if(p != r) {
      swap(rows[p], rows[r]);
      //swap(info[p], info[r]);
    }
    //info[r].second[r] = 1;
    for(int i = 0; i<n; i++){
      if(i != r && (rows[i]&bt(c))){
        rows[i] ^= rows[r];
        //info[i].second ^= info[r].second;
      }
    } r++;
  } return r;
}

B = [ n^(2/3), 1.26*n^(2/3)]
struct query{
  /* t = number of updates before this query*/
  int l, r, t, id;
  bool operator < (const query &x) const {
    if(l/B == x.l/B){
      if(r/B == x.r/B) return (t < x.t) ;
      return (r/B < x.r/B) ;
    }
  return l / B < x.l / B;
  }
};
struct upd{
  /// old = a[pos], a[pos] = cur
  int pos, old, cur;
};
void update(int pos, int x) {
  if (curL<=pos and pos<=curL)
    add(x), del(a[pos]);
  a[pos] = x;
}
t = totalupdates, curL = 1, curR = 0;
/// can start with t = 0 as well
for (int i = 1; i <= nq; i++) {
  int L = Q[i].l, R = Q[i].r, T = Q[i].t;
  while(t < T) t++, update(U[t].pos, U[t].cur);
  while(t > T) update(U[t].pos, U[t].old), t--;
  while(curL > L) add(a[--curL]);
  while(curR < R) add(a[++curR]);
  while(curL < L) del(a[curL++]);
  while(curR > R) del(a[curR--]);
  ans[Q[i].id] = something;
}

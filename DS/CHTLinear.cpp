// Minimum:
// M inc, x dec, useless(s-1, s-2, s-3)
// M dec, x inc, useless(s-3, s-2, s-1)
// Maximum:
// M inc, x inc, useless(s-3, s-2, s-1)
// M dec, x dec, useless(s-1, s-2, s-3)
//If queries are mot in order, use query2 O(logn).
typedef long long LL;
struct CHT {
  vector<LL> M; vector<LL> C; int ptr = 0;
  ///Use double comp if M,C is LL range
  bool useless(int l1, int l2, int l3) {
    return (C[l3]-C[l1])*(M[l1]-M[l2])
        <= (C[l2]-C[l1])*(M[l1]-M[l3]); }
  LL f(int id, LL x) { return M[id]*x+C[id]; }
  void add(LL m, LL c) {
    M.push_back(m); C.push_back(c);
    int s = M.size();
    while (s >= 3 && useless(s-3, s-2, s-1)) {
      M.erase(M.end()-2);C.erase(C.end()-2); s--;
    }
  }
  LL query(LL x) {
    if (ptr >= M.size()) ptr = M.size()-1;
    while (ptr < M.size()-1 && f(ptr, x)
           > f(ptr+1, x)) ptr++;/// > to < for max
    return f(ptr, x);
  }
  LL query2(LL x) {
    int lo=0, hi=M.size()-1;
    while(lo<hi) {
      int mid = (lo+hi)/2;
      /// change > to < for maximum
      if (f(mid, x) > f(mid+1, x))  lo = mid+1;
      else              hi = mid;
    } return f(lo, x);
  }
};

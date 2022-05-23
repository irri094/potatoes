// dp[n][k] = min(dp[i-1][k-1] + cost(i, n))
// using aliens trick, cost() is QF. O(n log^2)
namespace ArrayPartitionDP {
LL base_cost(int l, int r); ///define in code
long long C; int n;
PLL operator+ (const PLL &a, const PLL &b) {
  return PLL(a.first+b.first, a.second+b.second);
}
//Solves dp[i] = min(dp[j]+cost(j+1,i)), QF cost()
//returns {dp[n], min no of partitions}
PLL solve1D() {
  auto cost = [&](int l, int r) {
      return PLL(base_cost(l, r)+C, 1); };
  vector<PLL> dp(n+1); vector<int> opt(n+1);
  deque<pair<int, int>> dq; dq.push_back({0, 1});
  dp[0] = {0, 0};
  for (int i=1; i<=n; i++) {
    opt[i] = dq.front().first;
    dp[i] = dp[opt[i]] + cost(opt[i]+1, i);
    if (i == n) break;
    dq[0].second++;
    if (dq.size()>1 && dq[0].second==dq[1].second)
      dq.pop_front();
    int en = n;
    while(dq.size()) {
      int o=dq.back().first, st=dq.back().second;
      if (dp[o]+cost(o+1,st)>=dp[i]+cost(i+1,st))
        dq.pop_back();
      else {
        int lo = st, hi = en;
        while (lo < hi) {
          int mid = (lo+hi+1)/2;
          if (dp[o]+cost(o+1, mid) <
              dp[i]+cost(i+1, mid) )   lo = mid;
          else                         hi = mid-1;
        }
        if (lo < n) dq.push_back({i, lo+1});
        break;
      }
      en = st-1;
    }
    if (dq.empty()) dq.push_back({i, i+1});
  }
  return dp[n];
}
PLL check(LL c) { C = c; return solve1D();}
LL solve(int N, int k, LL lo, LL hi) {
  n = N;
  while (lo < hi) {
    LL mid = lo + (hi-lo)/2;
    if (check(mid).second > k)  lo = mid+1;
    else            hi = mid;
  }
  return check(lo).first - 1LL*k*lo;
} }

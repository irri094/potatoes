typedef long long LL;
struct PrimeSig{
  vector<int> primes;
  PrimeSig() {
    int MX = 100;
    vector<bool> isp(MX, 1);
    for (int i=2; i<MX; i++)
      if (isp[i]) {
        primes.push_back(i);
        for (int j=2*i; j<MX; j+=i) isp[j] = 0;
      }
  }
  LL LIM;
  vector<pair<vector<int>,LL>> ans; vector<int> ps;
  void go(int idx, LL val, int mx) {
    assert(ans.size() < 100000);
    assert(idx < primes.size());
    ans.push_back({ps, val});
    int p = primes[idx]; ps.push_back(0);
    for (int i=1; i<=mx; i++) {
      if (val > LIM/p)  break;
      ps.back()++; val *= p; go(idx+1, val, i);
    }
    ps.pop_back();
  }
  ///{signature, min value with signature} pair
  vector<pair<vector<int>, LL>> getAllSignature(LL lim) {
    LIM = lim; ans.clear(); ps.clear();
    go(0, 1, 100); return ans;
  }
};

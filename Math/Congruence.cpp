///Make mods const if possible
typedef long long LL; typedef pair<LL, LL> PLL;
LL power(LL a, LL b, LL m) {
  a = (a%m+m)%m; LL ans = 1;
  while (b) {
    if (b & 1) ans = (ans*a)%m;
    a = (a*a)%m; b >>= 1;
  }
  return ans;
}
LL egcd(LL a, LL b, LL &x, LL &y) {
  LL xx = y = 0; LL yy = x = 1;
  while (b) {
    LL q = a/b; LL t = b; b = a%b; a = t;
    t = xx; xx = x-q*xx; x = t;
    t = yy; yy = y-q*yy; y = t;
  }
  return a;
}
LL inverse(LL a, LL m) {
  LL x, y; LL g = egcd(a, m, x, y);
  if (g > 1) return -1; return (x%m+m)%m;
}
PLL CRT(LL m1, LL r1, LL m2, LL r2) {
  LL s, t; LL g = egcd(m1, m2, s, t);
  if (r1%g != r2%g) return PLL(0, -1);
  LL ss = ((s*r2)%m2)*m1, tt = ((t*r1)%m1)*m2;
  LL M = m1*m2, ans = ((ss+tt)%M+M)%M;
  return PLL(ans/g, M/g);
}
PLL CRT(const vector<LL> &m, const vector<LL> &r){
  PLL ans = PLL(r[0], m[0]);
  for (LL i = 1; i < m.size(); i++) {
    ans = CRT(ans.second, ans.first, m[i], r[i]);
    if (ans.second == -1) break;
  } return ans;
}
///computes x and y such that ax + by = c
bool LinearDiophantine(LL a,LL b,LL c,LL &x,LL &y){
  if (!a && !b) { if (c)      return false;
    x = y = 0;    return true;
  } if (!a) { if (c%b)    return false;
    x = 0; y = c/b; return true;
  } if (!b) { if (c%a)    return false;
    x = c/a; y = 0; return true;
  }
  LL g = gcd(a, b); if (c%g) return false;
  x = c/g * inverse(a/g, b/g); y = (c-a*x)/b;
  return true;
}
LL primitive_root(LL p) {
  if (p == 2) return 1;
  LL phi = p-1, n = phi; vector<LL> factor;
  for (int i=2; i*i<=n; ++i)
    if (n%i == 0) { factor.push_back(i);
      while (n%i==0)  n/=i;
    }
  if (n>1)  factor.push_back(n);
  for (int res=2; res<=p; ++res) {
    bool ok = true;
    for (int i=0; i<factor.size() && ok; ++i)
      ok &= power(res, phi/factor[i], p) != 1;
    if (ok)  return res;
  } return -1;
}
int discreteLog(int a, int b, int M) {
  a %= M, b %= M; int k = 1, add = 0, g;
  while ((g = gcd(a, M)) > 1) {
    if (b == k) return add; if (b % g)  return -1;
    b /= g, M /= g, ++add; k = (1LL*k*a/g)%M;
  }
  int RT = sqrt(M)+1, aRT = 1;
  for (int i=0; i<RT; i++) aRT = (aRT*1LL*a)%M;
  unordered_map<int, int> vals;
  for (int i=0, cur=b; i<=RT; i++) {
    vals[cur] = i; cur = (cur*1LL*a)%M;
  }
  for (int i=1, cur=k; i<=M/RT+1; i++) {
    cur = (cur*1LL*aRT)%M;
    if (vals.find(cur) != vals.end())
      return RT*i-vals[cur]+add;
  } return -1;
}
int discreteRoot(int a, int b, int P) {
  if (b%P == 0)   return  a == 0 ? -1 : 0;
  int g = primitive_root(P);
  int y = discreteLog(power(g, a, P), b, P);
  return y == -1 ? -1 : power(g, y, P);
}
vector<LL> berlekampMassey(vector<LL> s, int M) {
  if (s.empty()) return {};
  int n = s.size(), L = 0, m = 0;
  vector<LL> C(n), B(n), T;
  C[0] = B[0] = 1; LL b = 1;
  for (int i = 0; i < n; ++i) {
    ++m; LL d = s[i] % M;
    for (int j = 1; j <= L; ++j)
      d = (d + C[j] * s[i - j]) % M;
    if (!d) continue;
    T = C; LL coeff = d * power(b, M-2, M) % M;
    for (int j = m; j < n; ++j)
      C[j] = (C[j] - coeff * B[j - m]) % M;
    if (2*L > i) continue;
    L = i+1-L, B = T, b = d, m = 0;
  }
  C.resize(L + 1), C.erase(C.begin());
  for (LL &x : C) x = (M - x) % M;
  return C;
}

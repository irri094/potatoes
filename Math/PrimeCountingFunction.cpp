namespace PCF {//O(n^(3/4)/log). N = 10^13 -> 1.5s
const LL MAX = 1e13;
const int N = 7e6, M = 7, PM = 2*3*5*7*11*13*17;
bool isp[N];
int pr[N], pi[N]; int phi[M+1][PM+1], sz[M+1];
auto div = [](LL a,LL b)->LL{return double(a)/b;};
auto rt2 = [](LL x) -> int { return sqrtl(x); };
auto rt3 = [](LL x) -> int { return cbrtl(x); };
void init() {
  int cnt = 0; pi[0] = pi[1] = 0;
  for (int i=2; i<N; i++) isp[i] = true;
  for(int i=2; i<N; i++) {
    if(isp[i]) pr[++cnt] = i;
    pi[i] = cnt;
    for(int j=1; j<=cnt && i*pr[j]<N; j++) {
      isp[i * pr[j]] = false;
      if(i % pr[j] == 0)   break;
    }
  }
  sz[0] = 1;
  for(int i = 0; i <= PM; ++i)  phi[0][i] = i;
  for(int i = 1; i <= M; ++i) {
    sz[i] = pr[i] * sz[i - 1];
    for(int j = 1; j <= PM; ++j)
  phi[i][j]=phi[i-1][j]-phi[i-1][div(j,pr[i])];
  }
}
LL getphi(LL x, int s) {
  if(s == 0)  return x;
  if(s <= M)  return phi[s][x % sz[s]]
                    + (x/sz[s]) * phi[s][sz[s]];
  if(x <= 1LL*pr[s]*pr[s])   return pi[x] - s + 1;
  if(x <= 1LL*pr[s]*pr[s]*pr[s] && x < N) {
    int s2x = pi[rt2(x)];
    LL ans = pi[x]-(s2x+s-2)*(s2x-s+1)/2;
   for(int i=s+1;i<=s2x;++i)ans+=pi[div(x,pr[i])];
    return ans;
  }
  return getphi(x, s-1)-getphi(div(x,pr[s]), s-1);
}
LL getpi(LL x) {
  if(x < N)   return pi[x];
  LL ans = getphi(x, pi[rt3(x)]) + pi[rt3(x)] - 1;
  for(int i=pi[rt3(x)]+1,ed=pi[rt2(x)];i<=ed;++i)
    ans -= getpi(div(x,pr[i]))-i+1;
  return ans;
} }

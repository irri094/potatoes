namespace PSUM{ //for all powersum(k=1..N)
  const int N = 5e3 + 2;
  ll bern[N], sum[N], fac[N], ifac[N];
  ll bigmod(ll a, int n){
    ll ans = 1;
    while(n){
      if(n & 1)ans = ans * a % mod;
      n >>= 1; a = a * a % mod;}
    return ans;}
  void init(){//call me first
    fac[0] = 1;
    for(int i = 1; i < N; i++)
      fac[i] = fac[i - 1] * i % mod;
    ifac[N - 1] = bigmod(fac[N - 1], mod - 2);
    for(int i = N - 2; i >= 0; i--)
      ifac[i] = ifac[i + 1] * (i + 1) % mod;
    for(int i = 0; i < N; i++){ bern[i] = 1;
      for(int j = 0; j < i; j++){
        bern[i] = (bern[i]-fac[i] * ifac[j] % mod
         *ifac[i - j + 1] % mod * bern[j]) % mod;
        if(bern[i] < 0)bern[i] += mod;
      }}
  }
  //sum of i ^ k for 1 <= i <= n
  ll getPowerSum(ll n, int k){
    ll ans = 0, temp = n;
    for(int i=k;i>=0;i--,temp=temp*n%mod){
      ans = (ans + bern[i] * ifac[i] % mod * 
        ifac[k - i + 1] % mod * temp) % mod;
    }return (ll)ans * fac[k] % mod;
  }
  void build(ll n){
    init();
    for(int i = 0; i < N; i++)
      sum[i] = getPowerSum(n, i);}
}
// Find 1^k+2^k+...+n^k % mod 
/* x^k=sum(i=1 to k)Stirling2(k, i)* i!* ncr(x, i)
sum (x = 0 to n) x^k
=sum(i=0tok)Stirling2(k,i)*i!*sum(x=0ton)ncr(x,i)
=sum(i=0tok)Stirling2(k,i)*i!*ncr(n + 1, i + 1)
=sum(i=0tok)Stirling2(k,i)*i!*(n+1)!/(i+1)!/(n-i)!
=sum(i=0 to k)Stirling2(k,i)*(n - i + 1) *
		(n - i + 2) * ... (n + 1) / (i + 1)  */
ll S[105][105];
ll solve(int n, int k) {//(Shorter)
  S[0][0] = 1 % mod;
  for (int i = 1; i <= k; i++) {
    for (int j = 1; j <= i; j++) {
      if (i == j) S[i][j] = 1 % mod;
      else S[i][j] = (j * S[i - 1][j] +
	S[i - 1][j - 1]) % mod;}}
  ll ans = 0;
  for (int i = 0; i <= k; i++) {
    ll fact = 1, z = i + 1;
    for (ll j = n - i + 1; j <= n + 1; j++) {
			ll mul = j;
			if (mul % z == 0) {mul /= z;	z /= z;}
			fact = (fact * mul) % mod; 
		}
		ans = (ans + S[k][i] * fact) % mod;
	}
  return ans;
}

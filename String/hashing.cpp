/**Hash(ab..z) = a*p^n+b*p^(n-1)+... + z
Some Primes:
1000000007,1000000009,1000000861,1000099999(<2^30)
1088888881,1111211111,1500000001,1481481481(<2^31)
**/
typedef long long LL;typedef pair<LL, LL> PLL;
namespace Hashing {
  #define ff first
  #define ss second
  const PLL M = {1e9+7, 1e9+9};        
  const LL base = 1259; const int N = 1e6+7;                 
  PLL operator+ (const PLL& a, LL x)     
    {return {a.ff + x, a.ss + x};}
  PLL operator- (const PLL& a, LL x)     
    {return {a.ff - x, a.ss - x};}
  PLL operator* (const PLL& a, LL x)     
    {return {a.ff * x, a.ss * x};}
  PLL operator+ (const PLL& a, PLL x)    
    {return {a.ff + x.ff, a.ss + x.ss};}
  PLL operator- (const PLL& a, PLL x)    
    {return {a.ff - x.ff, a.ss - x.ss};}
  PLL operator* (const PLL& a, PLL x)    
    {return {a.ff * x.ff, a.ss * x.ss};}
  PLL operator% (const PLL& a, PLL m)    
    {return {a.ff % m.ff, a.ss % m.ss};}
  PLL pb[N];      ///powers of base mod M
  void hashPre() {  ///Call pre before everything
    pb[0] = {1,1};
    for(int i=1; i<N; i++)pb[i]=(pb[i-1]*base)%M;}
  ///Calculates hashes of all prefixes 
  vector<PLL> hashList(string s) {
    int n = s.size();
    vector<PLL> ans(n+1); ans[0] = {0,0};
    for (int i=1; i<=n; i++)    
      ans[i] = (ans[i-1] * base + s[i-1])%M;
    return ans;}
  ///Calculates hash of s[l..r](1 indexed)
  PLL subHash(const vector<PLL> &h,int l,int r)
  { return (h[r]+(M-h[l-1])*pb[r-l+1])%M; }
  ///Calculates Hash of a string
  PLL Hash (const string &s) {  PLL ans = {0,0};
      for (int i=0; i<s.size(); i++)  
        ans=(ans*base + s[i])%M;
      return ans;}
  PLL append(PLL cur, char c)
  {  return (cur*base + c)%M;  }
  PLL prepend(PLL cur, int k, char c)
  {  return (pb[k]*c + cur)%M;  }
  PLL replace(PLL cur, int i, char a, char b) 
  {  return cur + pb[i] * (M+b-a)%M;  }
  PLL pop_front(PLL hash, int len, char c) 
  {  return (hash + pb[len-1]*(M-c))%M;  }
  ///concatss two strings.k=length of the right
  PLL concat(PLL left, PLL right, int k) 
  {  return (left*pb[k] + right)%M;  }
  PLL power (const PLL& a, LL p) {
    PLL ans = {1, 1};
    while(p > 0){
      if (p&1)ans = (ans*a)%M;
      p >= 1; a = (a * a) % M;}
    return ans;}
  PLL inv(PLL a)  {
    if (M.ss == 1)  return power(a, M.ff-2);
    return power(a, (M.ff-1)*(M.ss-1)-1);}
  PLL invb = inv({base, base});
  PLL pop_back(PLL hash, char c) 
  {  return ((hash-c+M)*invb)%M;  }
  ///hash of str with size len repeated cnt times
  ///O(log n). For O(1), pre-calculate inverses
  PLL repeat(PLL hash, int len, LL cnt) {
    PLL mul=((pb[len*cnt]-1+M)*inv(pb[len]-1+M))%M;
    PLL ans = (hash*mul);
    if (pb[len].ff == 1)    ans.ff = hash.ff*cnt;
    if (pb[len].ss == 1)    ans.ss = hash.ss*cnt;
    return ans%M;}
}

typedef pair< int , int >PII;
ostream &operator<<(ostream &out, const PII &p) {
  return out<<"("<<p.first<<","<<p.second<<")";
}
const int maxn = 1e5+5; ///NOTICE
namespace DA {
  int wa[maxn],wb[maxn],wv[maxn],wc[maxn],r[maxn];
  int sa[maxn],rak[maxn],height[maxn],SIGMA=0;
  int cmp(int *r,int a,int b,int l) {
    return r[a] == r[b] && r[a+l] == r[b+l];
  }
  void da(int *r,int *sa,int n,int m) {
    int i,j,p,*x=wa,*y=wb,*t;
    for( i=0;i<m;i++) wc[i]=0;
    for( i=0;i<n;i++) wc[x[i]=r[i]]++;
    for( i=1;i<m;i++) wc[i] += wc[i-1];
    for( i= n-1;i>=0;i--)sa[--wc[x[i]]] = i;
    for( j= 1,p=1;p<n;j*=2,m=p){
      for(p=0,i=n-j;i<n;i++)y[p++] = i;
      for(i=0;i<n;i++)if(sa[i]>=j) y[p++]=sa[i]-j;
      for(i=0;i<n;i++)wv[i] = x[y[i]];
      for(i=0;i<m;i++) wc[i] = 0;
      for(i=0;i<n;i++) wc[wv[i]]++;
      for(i=1;i<m;i++) wc[i] += wc[i-1];
      for(i=n-1;i>=0;i--) sa[--wc[wv[i]]] = y[i];
      for(t=x,x=y,y=t,p=1,x[sa[0]]=0,i=1;i<n;i++)
        x[sa[i]]=cmp(y,sa[i-1],sa[i],j) ? p-1:p++;
    }
  }
  void calheight(int *r,int *sa,int n) {
    int i,j,k=0;
    for(i=1;i<=n;i++) rak[sa[i]] = i;
    for(i=0;i<n;height[rak[i++]] = k ) {
      for(k?k--:0,j=sa[rak[i]-1];
                              r[i+k]==r[j+k];k++);
    }
  }
  void suffixArray(const string &s,
      vector< int >&suffArray, vector< int >&lcp){
    int n = s.size(); SIGMA = 0;
    for(int i = 0; i < n; i ++) {
      if ('a'<=s[i]&&s[i]<='z') r[i] = s[i]-'a'+2;
      else r[i] = 1;  ///separators
      SIGMA = max(SIGMA, r[i]);
    }
    r[n] = 0; da(r, sa, n+1, SIGMA + 1);
    suffArray.resize(n);
    for(int i = 0; i<n; i++) suffArray[i]=sa[i+1];
    calheight(r,sa,n); lcp.resize(n-1);
    for (int i = 0; i+1<n; i++)lcp[i]=height[i+2];
  }
}
typedef vector<int>VI;const int K=20;int lg[maxn];
void pre() { ///CALL ME PLS
  lg[1] = 0;
  for (int i=2; i<maxn; i++) lg[i] = lg[i/2]+1;
}
struct RMQ{
  int N; VI v[K];
  RMQ(const VI &a) {
    N = a.size(); v[0] = a;
    for (int k = 0; (1<<(k+1)) <= N; k++) {
      v[k+1].resize(N);
      for (int i = 0; i-1+(1<<(k+1)) < N; i++) {
        v[k+1][i] = min(v[k][i], v[k][i+(1<<k)]);
      }
    }
  }
  int findMin(int i, int j) const {
    assert(i <= j); int k = lg[j-i+1];
    return min(v[k][i], v[k][j+1-(1<<k)]);
  }
};
PII extend(RMQ &rmq, int saSize, int ps, int len){
  int L = ps, R = ps;
  for (int k = K-1; k >= 0; k--) {
    int r = R+(1<<k); if (r >= saSize) continue;
    if (rmq.findMin(ps, r-1) >= len) R = r;
  }
  for (int k = K-1; k >= 0; k--) {
    int l = L-(1<<k); if (l < 0) continue;
    if (rmq.findMin(l, ps-1) >= len) L = l;
  }
  return PII(L, R);
}
///len retrnd by backstep() must b min-ed from out
struct BackStepper {
  vector< int >startsWith[26];
BackStepper(const string &s,const vector<int>&sa){
    for (int i = 0; i < sa.size(); i++) {
      if (sa[i] > 0)
        startsWith[s[sa[i]-1]-'a'].push_back(i);
    }
    startsWith[s.back()-'a'].push_back(s.size());
  }
  /**Return <len, j> s.t. s[j] = c and suffix[j+1]
     shares the longest prefix with suffix[i]
      Returns <0, 0> if no such index exists. */
  PII backstep(int i,int c,const vector< int >&sa,
          const vector< int >&ra, const RMQ &rmq){
    if (startsWith[c].empty()) return PII(0, 0);
    int ri = ra[i];
    int idx = lower_bound(startsWith[c].begin(),
    startsWith[c].end(),ri)-startsWith[c].begin();
    if (idx < startsWith[c].size() &&
      startsWith[c][idx]==ri) { ///same pos again
      return PII(ra.size()-i+1, i-1);
    }
    PII rt(-1, -1);
    if (idx > 0) {
      int ci = startsWith[c][idx-1];
      rt = PII(rmq.findMin(ci, ri-1)+1, sa[ci]-1);
    }
    if (idx < startsWith[c].size()) {
      int ci = startsWith[c][idx];
      if (ci==sa.size())
        rt=max(rt,PII(1, sa.size()-1));
      else rt = max(rt,
        PII(rmq.findMin(ri, ci-1)+1, sa[ci]-1));
    }
    return rt;
  }
};
/** n -> string length, SZ -> number of nodes in
suff tree, 1 is root (the empty string), length of
node u (starting from root) is length[u], for an
edge between node u and it's child v, edge length
is length[v]-length[u], node u belongs to all of
[L[u], R[u]] suffixes in Suffix Array */
const int MXND = maxn*2+7;
struct SuffixTree {
  int nxt[maxn]; /// nxt[i] = position of next #
  vector<int>edg[MXND], leaves[MXND];
  int lnth[MXND],L[MXND],R[MXND], SZ, n;
  void buildGraph(const vector<int>&sa,
      vector<int>&lcp) {
    //edg and leaves must be cleared if test cases
    SZ = 0; vector<int> stk{++SZ};
    n = sa.size(), L[SZ] = 0, R[SZ] = n-1;
    lnth[SZ] = 0;lcp.push_back(0);int last = -1;
    for(int i = 0, sf = 1; i+sf<=n; i+=sf, sf^=1){
      int left=i-(sf^1);//sf=suflen/lcp being used
      //int curlcp=(sf)?n-sa[i]:(i?lcp[i-1]:0);
      /// ^^ single string;
      int curlcp = (sf)?nxt[sa[i]]-sa[i]:
        (i?lcp[i-1]:0);  /// multiple string
      while(curlcp < lnth[stk.back()]){
        R[stk.back()]=i-(sf^1),left=L[stk.back()];
        last = stk.back(), stk.pop_back();
        if(curlcp <= lnth[stk.back()])
          edg[stk.back()].push_back(last),last=-1;
      }
      if(curlcp > lnth[stk.back()]){
        stk.push_back(++SZ);
      if(last!=-1)edg[SZ].push_back(last),last=-1;
        lnth[SZ]=curlcp,L[SZ]=left;
      }
    }
  }
  void buildLeaves(const vector< int >&sa){
    for(int i = 1; i<=SZ; i++){
      int r=(edg[i].empty()?R[i]:L[edg[i][0]]-1);
      for(int j = L[i]; j<=r; j++){
        leaves[i].push_back(sa[j]);
      }
    }
  }
} st;
void dfs(const string &s,vector<int>&sa,int u){
  for (int v : st.edg[u]) {
    cout << u << " -:";
    cout << s.substr(sa[st.L[v]]+st.lnth[u],
        st.lnth[v]-st.lnth[u]);
    cout << ":- " << v; cout << " {";
    for (int x : st.leaves[v]) cout << " " << x;
    cout << " }"; cout << endl; dfs(s, sa, v);
  } ///suffLen[i] denotes the actual
} /// suffix length of sa[i]
int suffLen[maxn];
void sanitize(const string &s,
    const vector<int>&sa, vector< int >&lcp) {
  int n = s.size(); st.nxt[n] = n;
  for (int i = n-1; i >= 0; i--) {
   st.nxt[i]=st.nxt[i+1];if(s[i]=='#')st.nxt[i]=i;
  }
  for (int i = 0; i < n; i++) {
    suffLen[i] = st.nxt[ sa[i] ] - sa[i];
  }
  for (int i = 0; i+1 < n; i++) {
  lcp[i]=min(lcp[i],min(suffLen[i],suffLen[i+1]));
  }
}
int main() {
  pre();///MUST for RMQ to be working
  string s = "abac#ababa";
  vector<int>sa,lcp; DA::suffixArray(s, sa, lcp);
  sanitize(s,sa,lcp); cout<<"suffix order:"<<endl;
  for(int i:sa)cout << s.substr(i) << endl;
  cout << "lcp:"; for(int i:lcp) cout << " " << i;
  cout << endl;st.buildGraph(sa, lcp);
  st.buildLeaves(sa);dfs(s, sa, 1);return 0;
}

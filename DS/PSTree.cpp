#include<bits/stdc++.h>
using namespace std;
typedef pair<int,int>PII;
typedef long long LL;
#define all(x) x.begin(), x.end()
const int MAXN = 2e5+7;
const int LOGN = 19; ///ATTENTION
namespace PSTree {
  struct PSNode { int cnt, lc, rc; };
  PSNode tr[MAXN*LOGN+100];
  ///since there is no build function, you either
  ///write an appropriate build function or make
  ///sure cnt = identity; zero-th node as initial
  int counter;
  void clear() { counter = 0; }
  int update(int u, int l, int r, int idx, int v){
    if (idx < l || r < idx) return u;
    if (l == r) {
      int nd = ++counter; tr[nd].cnt=tr[u].cnt+v;
      return nd;
    }
    int mid = (l+r)/2; int nd = ++counter;
    tr[nd].lc =update(tr[u].lc, l, mid, idx, v);
    tr[nd].rc =update(tr[u].rc, mid+1, r, idx, v);
   tr[nd].cnt=tr[tr[nd].lc].cnt+tr[tr[nd].rc].cnt;
    return nd;
  }
  int getKth(vector<PII>vp, int l, int r, int k){
    if (l==r) return l; int mid = (l+r)/2,bam = 0;
    for (PII pr : vp) {
      int lc = tr[pr.first].lc;
      bam += tr[lc].cnt * pr.second;
    }
    if (k <= bam) {
      for (PII &pr:vp) pr.first = tr[pr.first].lc;
      return getKth(vp, l, mid, k);
    } else {
      for (PII &pr:vp) pr.first = tr[pr.first].rc;
      return getKth(vp, mid+1, r, k-bam);
    }
  }
}
int ar[MAXN], rt[MAXN];
int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);
  int n, m;
  cin >> n >> m;
  vector< int >vas(n);
  for (int i = 1; i <= n; i++) {
    cin >> ar[i];
    vas[i-1] = ar[i];
  }
  sort(all(vas));
  vas.erase(unique(all(vas)), vas.end());
  int sz = vas.size();
  for (int i = 1; i <= n; i++) {
    ar[i]=lower_bound(all(vas),ar[i])-vas.begin();
    rt[i]=PSTree::update(rt[i-1],0,sz-1,ar[i],1);
  }
  while (m--) {
    int l, r, k;
    cin >> l >> r >> k;
    int v = PSTree::getKth({PII(rt[r],+1),
                    PII(rt[l],-1)}, 0, sz-1, k+1);
    cout << vas[v] << "\n";
  }
  return 0;
}

#include<bits/stdc++.h>
using namespace std;
const int MAXN = 1e5+7; /// change this
struct nd { int mn, cnt; };
nd tr[4*MAXN]; int lazy[4*MAXN];
///1. Merge left and right
nd combine(const nd &a, const nd &b) {
  if (a.mn < b.mn) return a;
  else if (a.mn > b.mn) return b;
  nd rt;rt.mn=a.mn;rt.cnt=a.cnt+b.cnt;return rt;
}
///2. Push lazy down and merge lazy
void propagate(int u, int l, int r) {
  if (lazy[u]) {
    tr[u].mn += lazy[u];
    if (l != r) {
      lazy[u*2] += lazy[u];lazy[u*2+1] += lazy[u];
    }
    lazy[u] = 0;
  }
}
int a[MAXN];
void build(int u, int l, int r) {
  lazy[u] = 0;     ///3. Initialize
  if (l==r) {
    tr[u].mn = 0; tr[u].cnt = a[l]; return;
  }
  int mid = (l+r)/2;
  build(u*2, l, mid); build(u*2+1, mid+1, r);
  tr[u] = combine(tr[u*2], tr[u*2+1]);
}
void update(int u,int l,int r,int x,int y,int v){
  propagate(u, l, r);
  if (r < x || y < l) return;
  if (x <= l && r <= y) {
    lazy[u] += v; ///4. Merge lazy
    propagate(u, l, r); return;
  }
  int mid = (l+r)/2; update(u*2,l,mid,x,y,v);
  update(u*2+1,mid+1,r,x,y,v);
  tr[u] = combine(tr[u*2], tr[u*2+1]);
}
nd query(int u, int l, int r, int x, int y) {
  propagate(u, l, r);
  if (x <= l && r <= y) return tr[u];
  int mid = (l+r)/2;
  if (y <= mid) return query(u*2, l, mid, x, y);
  if (mid < x)return query(u*2+1, mid+1, r, x, y);
  return combine(query(u*2, l, mid, x, y),
                 query(u*2+1, mid+1, r, x, y));
}

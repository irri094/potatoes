#include<bits/stdc++.h>
using namespace std;
typedef long long LL;
mt19937 rng(chrono::steady_clock::now().
            time_since_epoch().count());
typedef struct item * pitem;
struct item {
  int prior, value, cnt; LL sum; bool rev;
  item(int value):prior(rng()), value(value) {
    cnt = 0;rev = 0;sum = value;l = r = nullptr;
  }
  pitem l, r;
};
/// 0-indexed
namespace Treap {
  int cnt(pitem it){return it!=nullptr?it->cnt:0;}
  LL sum(pitem it) {return it!=nullptr?it->sum:0;}
  void upd_cnt (pitem it) {
    if (it!=nullptr) {
      it->cnt=cnt(it->l)+cnt(it->r)+1;
      it->sum=sum(it->l)+sum(it->r)+it->value;
    }
  }
  void push (pitem it) {
    if (it!=nullptr && it->rev==true){
      it->rev = false; swap (it->l, it->r);
      if (it->l)  it->l->rev ^= true;
      if (it->r)  it->r->rev ^= true;
    }
  }
  void merge (pitem & t, pitem l, pitem r) {
    push (l); push (r);
    if (l==nullptr || r==nullptr)
      t = (l!=nullptr) ? l : r;
    else if (l->prior > r->prior)
      merge (l->r, l->r, r),  t = l;
    else merge (r->l, l, r->l),  t = r;
    upd_cnt (t);
  }
  /// r will have key-th and above
  void split(pitem t,pitem &l,pitem &r,int key,
      int add = 0) {
    if (t==nullptr) { l = r = nullptr; return; }
    push(t); int cur_key = add + cnt(t->l);
    if (key <= cur_key)
      split (t->l, l, t->l, key, add),  r = t;
    else
      split(t->r,t->r,r,key,add+1+cnt(t->l)), l=t;
    upd_cnt (t);
  }
  void reverse (pitem &t, int l, int r) {
    pitem t1, t2, t3; split (t, t1, t2, l);
    split (t2, t2, t3, r-l+1); assert(t2 != NULL);
    t2->rev^=true; merge(t,t1,t2); merge(t,t,t3);
  }
  LL query(pitem &t, int l, int r) {
    pitem t1, t2, t3; split(t, t1, t2, l);
    split(t2, t2, t3, r-l+1); LL ans = t2->sum;
    merge(t,t1,t2); merge(t,t,t3); return ans;
  }
  /// before key'th element
  void insert (pitem & t, int key, int value) {
    pitem x = new item(value); pitem L, R;
    split(t, L, R, key); merge(L, L, x);
    merge(t, L, R); upd_cnt(t);
  }
  int erase (pitem & t, int key) {
    assert(cnt(t) > key); pitem L, MID, R;
    split(t, L, MID, key); split(MID, MID, R, 1);
    merge(t, L, R); upd_cnt(t);
    int rt = MID->value; delete MID; return rt;
  }
  void output (pitem t, vector< int >&v) {
    if (t==nullptr) return;
    push (t); output (t->l, v);
    v.push_back(t->value); output (t->r, v);
  }
}
int main() {
  int n, q, m; cin >> n >> q >> m;
  pitem tr = nullptr;
  for (int i = 0; i < n; i++) {
    int x; cin >> x;
    Treap::insert(tr, i, x);
  }
  while (q--) {
    int t, l, r;
    cin >> t >> l >> r; l--; r--;
    if (t==1) {
      int x = Treap::erase(tr, r);
      Treap::insert(tr, l, x);
    } else { Treap::reverse(tr, l, r); }
  }
  vector< int >v; Treap::output(tr, v);
  while (m--) {
    int i; cin >> i;
    cout << v[i-1] << " ";
  }
}



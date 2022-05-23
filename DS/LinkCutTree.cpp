#include<bits/stdc++.h>
using namespace std;
const int MOD = 998244353;
int sum(int a, int b) {
  return a+b >= MOD ? a+b-MOD : a+b;
}
int mul(int a, int b) {
  return (a*1LL*b)%MOD;
}
typedef pair< int , int >Linear;
Linear compose(const Linear &p, const Linear &q) {
  return Linear(mul(p.first, q.first),
    sum(mul(q.second, p.first), p.second));
}
struct SplayTree {
  struct Node {
    int ch[2] = {0, 0}, p = 0;
    long long self = 0, path = 0;//Path aggregates
    long long sub = 0, vir = 0;//Subtree aggregate
    int size = 1; bool flip = 0;// Lazy tags
    Linear _self{1, 0}, shoja{1, 0}, ulta{1, 0};
  };
  vector<Node> T;
  SplayTree(int n) : T(n + 1) { T[0].size = 0; }
  void push(int x) {
    if (!x || !T[x].flip) return;
    int l = T[x].ch[0], r = T[x].ch[1];
    T[l].flip ^= 1, T[r].flip ^= 1;
    swap(T[x].ch[0], T[x].ch[1]); T[x].flip = 0;
    swap(T[x].shoja, T[x].ulta);
  }
  void pull(int x) {
    int l=T[x].ch[0],r=T[x].ch[1];push(l);push(r);
    T[x].size = T[l].size + T[r].size + 1;
    T[x].path = T[l].path + T[x].self + T[r].path;
    T[x].sub=T[x].vir+T[l].sub+T[r].sub+T[x].self;
    T[x].shoja = compose(T[r].shoja,
                compose(T[x]._self, T[l].shoja));
    T[x].ulta  = compose(T[l].ulta,
                compose(T[x]._self, T[r].ulta));
  }
  void set(int x, int d, int y) {
    T[x].ch[d] = y; T[y].p = x; pull(x);
  }
  void splay(int x) {
    auto dir = [&](int x) {
      int p = T[x].p; if (!p) return -1;
      return T[p].ch[0]==x?0:T[p].ch[1]==x?1:-1;
    };
    auto rotate = [&](int x) {
      int y = T[x].p,z=T[y].p,dx=dir(x),dy=dir(y);
      set(y, dx, T[x].ch[!dx]); set(x, !dx, y);
      if (~dy) set(z, dy, x); T[x].p = z;
    };
    for (push(x); ~dir(x); ) {
      int y = T[x].p,z = T[y].p; push(z); push(y);
      push(x); int dx = dir(x), dy = dir(y);
      if (~dy) rotate(dx!=dy?x:y);  rotate(x);
    }
  }
  int KthNext(int x, int k) {
    assert(k > 0); splay(x); x = T[x].ch[1];
    if (T[x].size < k) return -1;
    while (true) {
      push(x); int l = T[x].ch[0], r = T[x].ch[1];
      if (T[l].size+1 == k) return x;
      if (k <= T[l].size) x = l;
      else k -= T[l].size+1, x = r;
    }
  }
};
struct LinkCut : SplayTree {
  LinkCut(int n) : SplayTree(n) {}
  int access(int x) {
    int u = x, v = 0;
    for (; u; v = u, u = T[u].p) {
      splay(u); int& ov = T[u].ch[1];
      T[u].vir += T[ov].sub; T[u].vir -= T[v].sub;
      ov = v; pull(u);
    }
    splay(x); return v;
  }
  void reroot(int x) {
    access(x); T[x].flip ^= 1; push(x);
  }
  ///makes v parent of u !(u must be a root)
  void Link(int u, int v) {
    reroot(u); access(v); T[v].vir += T[u].sub;
    T[u].p = v; pull(v);
  }
  ///removes edge between u and v
  void Cut(int u, int v) {
    int _u = FindRoot(u); reroot(u); access(v);
    T[v].ch[0] = T[u].p = 0; pull(v); reroot(_u);
  }
  //Rooted tree LCA.Returns 0 if u v not connected
  int LCA(int u, int v) {
    if (u == v) return u; access(u);
    int ret = access(v); return T[u].p ? ret : 0;
  }
  //Query subtree of u where v is outside the sbtr
  long long Subtree(int u, int v) {
    int _v = FindRoot(v); reroot(v); access(u);
    long long ans = T[u].vir + T[u].self;
    reroot(_v); return ans;
  }
  long long Path(int u, int v) {
    int _u = FindRoot(u); reroot(u); access(v);
    long long ans = T[v].path; reroot(_u);
    return ans;
  }
  Linear _Path(int u, int v) {
    reroot(u); access(v); return T[v].shoja;
  }
  void Update(int u, long long v) {
    access(u); T[u].self = v; pull(u);
  }
  void _Update(int u, Linear v) {
    access(u); T[u]._self = v; pull(u);
  }
  int FindRoot(int u) {
    access(u);
    while (T[u].ch[0]) { u = T[u].ch[0]; push(u);}
    access(u); return u;
  }
  ///k-th node (0-indexed) on the path from u to v
  int KthOnPath(int u, int v, int k) {
    if (u == v) return k == 0 ? u : -1;
    int _u = FindRoot(u); reroot(u); access(v);
    int ans = KthNext(u, k); reroot(_u);
    return ans;
  }
};
int main() {
  int n, q; cin >> n >> q; LinkCut lct(n);
  for (int i = 1; i <= n; i++) {
    Linear l; cin >> l.first >> l.second;
    lct._Update(i, l);
  }
  for (int i = 1; i < n; i++) {
    int u, v; cin >> u >> v; lct.Link(u+1, v+1);
  }
  while (q--) {
    int op; cin >> op;
    if (op == 0) {
      int u, v, w, x; cin >> u >> v >> w >> x;
      lct.Cut(u+1, v+1); lct.Link(w+1, x+1);
    } else if (op == 1) {
      int p; Linear l; cin>>p>>l.first>>l.second;
      lct._Update(p+1, l);
    } else {
      int u, v, x; cin >> u >> v >> x;
      Linear l = lct._Path(u+1, v+1);
      cout<<sum(mul(l.first, x), l.second)<<"\n";
    }
  }
}

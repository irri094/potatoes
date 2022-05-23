/// Min value query
struct func{
  ll operator()(ll x){...}
} tree[4*nmax], maxfunc;
#define lc (id<<1)
#define rc ((id<<1)|1)
void build(int id, int l, int r){
  tree[id] = maxfunc; if(l+1 == r) return;
  int mid = (l+r)/2;
  build(lc, l, mid); build(rc, mid, r);
}
void add_func(int id, int l, int r, func f){
  int mid = (l+r)/2;
  bool lefbad = f(l) < tree[id](l);
  bool midbad = f(mid) < tree[id](mid);
  if(midbad) swap(f, tree[id]);
  if(l + 1 == r) return;
  else if(lefbad!=midbad) add_func(lc, l, mid, f);
  else add_func(rc, mid, r, f);
}
ll get_val(int id, int l, int r, ll x){
  ll tmp = tree[id](x);
  if(l + 1 == r) return tmp;
  int mid = (l+r)/2;
  if(x < mid) return min(get_val(lc,l,mid,x),tmp);
  else return min(get_val(rc, mid, r, x), tmp);
}

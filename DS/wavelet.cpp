int a[N];
struct weblet{
  int lo, hi;
  weblet *l=0, *r=0;
  vector<int> b, c;
// call weblet(a+1, a+n+1, minval, maxval)
  weblet(int *from, int *to, int x, int y){
    lo = x, hi = y;
    if( from >= to) return;
    if( hi == lo ){
      b.reserve(to-from+1), b.pb(0);
      c.reserve(to-from+1), c.pb(0);
      for(auto it=from; it!=to; it++){
        b.pb(b.back() + 1);
        c.pb(c.back()+*it);
      }
      return;
    }
    int mid = (lo+hi)/2;
    auto f = [mid](int x){
        return x <= mid;
    };
    b.reserve(to-from+1), b.pb(0);
    c.reserve(to-from+1), c.pb(0);
    for(auto it = from; it != to; it++){
      b.pb(b.back() + f(*it));
      c.pb(c.back() + *it);
    }
    auto pivot = stable_partition(from, to, f);
    l = new weblet(from, pivot, lo, mid);
    r = new weblet(pivot, to, mid+1, hi);
  }
  void swapadjacen(int i){/// i with i+1{
    if(lo == hi) return ;
    b[i]= b[i-1] + b[i+1] - b[i];
    c[i] = c[i-1] + c[i+1] - c[i];
    if( b[i+1]-b[i] == b[i] - b[i-1]){
      if(b[i] -b[i-1])
        return this->l->swapadjacen(b[i]);
      else return this->r->swapadjacen(i-b[i]);
    }
    else return ;
  }
  int kth(int l, int r, int k){
    if(l > r) return 0;
    if(lo == hi) return lo;
    int inLeft = b[r] - b[l-1];
    int lb = b[l-1], rb = b[r];
    if(k<=inLeft)
      return this->l->kth(lb+1, rb, k);
    return this->r->kth(l-lb, r-rb, k-inLeft);
  }
  int LTE(int l, int r, int k){
    if(l > r or k < lo) return 0;
    if(hi <= k) return r - l + 1;
    int lb = b[l-1], rb = b[r];
    return this->l->LTE(lb+1, rb, k)+
      this->r->LTE(l-lb, r-rb, k);
  }
  int count(int l, int r, int k){
    if(l > r or k < lo or k > hi) return 0;
    if(lo == hi) return r - l + 1;
    int lb = b[l-1], rb = b[r], md = (lo+hi)/2;
    if(k<=md) return this->l->count(lb+1,rb,k);
    return this->r->count(l-lb, r-rb, k);
  }
  int sumk(int l, int r, int k){/// sumof <=k
    if(l > r or k < lo) return 0;
    if(hi <= k) return c[r] - c[l-1];
    int lb = b[l-1], rb = b[r];
    return this->l->sumk(lb+1, rb, k) +
      this->r->sumk(l-lb, r-rb, k);
  }
  ~weblet(){
    if(l) delete l;
    if(r) delete r;
  }
};

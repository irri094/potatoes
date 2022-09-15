#ifdef LOCAL
#define Gene template< class
#define Rics printer& operator,
Gene c> struct rge{c b, e;};
Gene c> rge<c> range(c i, c j){ return {i, j};}
struct printer{
  ~printer(){cerr<<endl;}
  Gene c >Rics(c x){ cerr<<boolalpha<<x; return
  *this;}
  Rics(string x){cerr<<x;return *this;}
  Gene c, class d >Rics(pair<c, d> x){ return
  *this,"(",x.first,", ",x.second,")";}
  Gene ... d, Gene ...> class c >Rics(c<d...> x)
  { return *this, range(begin(x), end(x));}
  Gene c >Rics(rge<c> x){
  *this,"["; for(auto it = x.b; it != x.e; ++it)
  *this,(it==x.b?"":", "),*it; return *this,"]";}
};
#define debug() cerr<<"LINE "<<__LINE__\
<<" >> ", printer()
#define dbg(x) "[",#x,": ",(x),"] "
#define tham getchar()
#endif
mt19937_64 rng((unsigned)chrono::system_clock\
::now().time_since_epoch().count());
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt")

#include<bits/stdc++.h>
using namespace std;
typedef long long LL;
template< typename T >
struct SlopeTrick {
  T INF = numeric_limits< T >::max() / 3;
  T min_f, add_l, add_r;
  priority_queue< T, vector< T >, less<> > L;
  priority_queue< T, vector< T >, greater<> > R;
  void push_R(const T &a) { R.push(a - add_r); }
  T top_R() const { return R.top() + add_r; }
  T pop_R() { T val=top_R(); R.pop(); return val;}
  void push_L(const T &a) { L.push(a - add_l); }
  T top_L() const { return L.top() + add_l; }
  T pop_L() { T val=top_L(); L.pop(); return val;}
public:
  SlopeTrick() : min_f(0), add_l(0), add_r(0) {
    L.push(-INF); R.push(INF);
  }
  // f(x) += a
  void add_all(const T &a) {
    min_f += a;
  }
  // add \_ ; f(x) += max(a - x, 0)
  void add_a_minus_x(const T &a) {
    if (a > top_R()) {
      min_f+=a-top_R(); push_L(pop_R());push_R(a);
    } else {
      push_L(a);
    }
  }
  // add _/ ; f(x) += max(x - a, 0)
  void add_x_minus_a(const T &a) {
    if (top_L() > a) {
      min_f+=top_L()-a;push_R(pop_L());push_L(a);
    } else {
      push_R(a);
    }
  }
  // add \/ ; f(x) += max(x - a, 0)
  void add_abs(const T &a) {
    add_a_minus_x(a); add_x_minus_a(a);
  }
  // \/ -> \_ ; f_{new} (x) = min f(y) (y <= x)
  void clear_right(){ while(R.size()>=2)R.pop(); }
  // \/ -> _/ ; f_{new} (x) = min f(y) (y >= x)
  void clear_left(){ while(L.size()>=2) L.pop(); }
  // \/. -> .\/ ; f_{new} (x) = f(x - a)
  void shift(const T &a) { add_l+=a; add_r+=a; }
  T get(const T &x) {
    T ret = min_f;
    if (!L.empty() && x < top_L()) {
      while(!L.empty())ret += max(T(0),pop_L()-x);
    }
    if (!R.empty() && top_R() < x) {
      while(!R.empty())ret += max(T(0),x-pop_R());
    }
    return ret;
  }
};
void SmallToLarge(SlopeTrick<LL> &from,
                  SlopeTrick<LL>&to) {
  if (from.L.size()+from.R.size()>
      to.L.size()+to.R.size()) swap(from, to);
  while (from.L.size() >= 2) {
    to.add_a_minus_x(from.pop_L());
  }
  while (from.R.size() >= 2) {
    to.add_x_minus_a(from.pop_R());
  }
  to.min_f += from.min_f;
}
const int MAXN = 3e5+7; int P[MAXN], C[MAXN];
vector<int>child[MAXN];SlopeTrick<LL> trick[MAXN];
int main() {
  ios::sync_with_stdio(false);
  cin.tie(0);
  int n, m; cin >> n >> m;
  for (int i = 2; i <= n+m; i++) {
    cin >> P[i] >> C[i];child[P[i]].push_back(i);
  }
  for (int i = n+1; i <= n+m; i++) {
    trick[i].add_abs(C[i]);
  }
  for (int i = n; i > 0; i--) {
    for (int c : child[i]) {
      SmallToLarge(trick[c], trick[i]);
    }
    { ///clearing all slopes greater than 1
      LL save = trick[i].top_R();
      trick[i].clear_right();
      trick[i].push_R(save);
    }
    if (i > 1) {
      trick[i].push_L(trick[i].pop_L()+C[i]);
      trick[i].push_R(trick[i].pop_R()+C[i]);
    }
  }
  cout << trick[1].min_f << "\n";
  return 0;
}

///0-based indexing
///p[0][i] = length of half palin around hlf ind i
///p[1][i] = length of half palin around char i
struct Manacher {
  vector<int> p[2];
  Manacher(const string &s) {
    int n = s.size();
    p[0].resize(n+1); p[1].resize(n);
    for (int z=0; z<2; z++)
      for (int i=0, l=0, r=0; i<n; i++) {
        int t = r-i+!z;
        if (i<r) p[z][i] = min(t, p[z][l+t]);
        int L = i-p[z][i], R = i+p[z][i]-!z;
        while (L>=1 && R+1<n && s[L-1] == s[R+1])
          p[z][i]++, L--, R++;
        if (R>r) l=L, r=R;
      }
  }
  bool ispalin(int l, int r) {
    int mid = (l+r+1)/2, sz = r-l+1;
    return 2*p[sz%2][mid]+sz%2 >= sz;
  }
};

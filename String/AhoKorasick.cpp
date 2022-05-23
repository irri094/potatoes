const int MAXX = 1e6+7;  //total length of strings
namespace AhoCorasick {
  const int sigma = 26, offset = 'a';
  struct Vertex {
    int next[sigma];    /// indices of child node
    int lvl = 0;        /// depth of the node
    bool leaf = false;  /// if it is a last char
    int p = -1;         /// index of parent node
    char pch;           /// parent character
    int link = -1;      /// suffix link for vertex
    int go[sigma];      /// where to go from here
    Vertex(int p=-1, char ch='$'): p(p), pch(ch) {
      fill(next, next+sigma, -1);
      fill(go, go+sigma, -1);
    }
  } t[MAXX];
  int exit_link[MAXX]; int used = 1;
  void new_test_case() {
      used = 1; t[0] = Vertex();
  }
  int add_string(string const& s) {
    int v = 0;
    for (char ch : s) {
      int c = ch - offset;
      if (t[v].next[c] == -1) {
        t[v].next[c]=used; t[used]=Vertex(v,ch);
        t[used].lvl = t[v].lvl+1; used++;
      }
      v = t[v].next[c];
    }
    t[v].leaf = true; return v;
  }
  void build() {
    queue< int >q; q.push(0);
    while (!q.empty()) {
      int v = q.front(); q.pop();
      {   ///calculate suffix link
        if (v == 0 || t[v].p == 0) t[v].link = 0;
        else t[v].link =
          t[t[t[v].p].link].go[t[v].pch-offset];
      }
      {   ///calculate exit link
        if (v == 0 || t[v].p == 0) exit_link[v]=0;
        else if (t[t[v].link].leaf)
          exit_link[v] = t[v].link;
        else exit_link[v]=exit_link[t[v].link];
      }
      for (int i = 0; i < sigma; i++) {
        if (t[v].next[i]!=-1)q.push(t[v].next[i]);
      }
      /// save go values in next to save memory!
      for (int c = 0; c < sigma; c++) {
        if (t[v].next[c] != -1)
          t[v].go[c] = t[v].next[c];
        else
          t[v].go[c] = v==0? 0:t[t[v].link].go[c];
      }
    }
  }
}
int main() {
    int n, k; cin >> n >> k;
    AhoCorasick::new_test_case();
    for (int i = 1; i <= k; i++) {
      string t;cin>>t; AhoCorasick::add_string(t);
    }
    AhoCorasick::build();
}

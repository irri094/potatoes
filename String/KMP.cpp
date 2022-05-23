const int ALPHA = 26;
/// builds the prefix automaton in O(N*ALPHA)
vector< vector< int > >automaton;
void buildAutomaton(const string& s) {
  automaton.clear(); int n = s.size(), k = 0;
  for (int i = 0; i <= n; i++)
    automaton.emplace_back(ALPHA, 0);
  automaton[0][s[0]-'a'] = 1;
  for (int i = 1; i <= n; i++) {
    automaton[i] = automaton[k];
    if (i < n) {
      automaton[i][s[i]-'a'] = i+1;
      k = automaton[k][s[i]-'a'];
    }
  }
}
vector< int >prefixFunction(const string& s) {
  int n = s.size(), k = 0; /// 1-indexed
  vector< int >v(n+1); v[1] = 0;
  for (int i = 2; i <= n; i++) {
    while (k > 0 && s[k] != s[i-1]) k = v[k];
    if (s[k]==s[i-1]) k++;
    v[i] = k;
  }
  return v;
}
int matcher(const string& txt,const string& ptrn){
  vector< int >pi = prefixFunction(ptrn);
  int matchCount = 0, k = 0;
  for (int i = 0; i < txt.size(); i++) {
    while (k > 0 && txt[i]!=ptrn[k]) k = pi[k];
    if (txt[i]==ptrn[k]) k++;
    if (k==ptrn.size()) {
      matchCount++; k = pi[k];
    }
  }
  return matchCount;
}

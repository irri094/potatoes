const int MAXN = 1e5+7; ///max length of string+3
namespace PalTree {
  struct node {
    int len; ///length of the pal of this node
    int sufflink; ///largest suff pal of this node
    int chain;///#of nodes on chain of suff links
    int next[26];///next[c] is the pal by adding c
  } tr[MAXN];
  int size;///# of nodes currently in Pal tr
  int suff;///max suff pal of curr prcessed prefix
  string s;///string we will built our Paltr on
  bool addLetter(int pos) {
    int cur = suff, curlen = 0, let = s[pos]-'a';
    while (true) {
      curlen = tr[cur].len;
      if(pos-1-curlen>=0&&s[pos-1-curlen]==s[pos])
        break; cur = tr[cur].sufflink;
    }
    if (tr[cur].next[let]) {
      suff = tr[cur].next[let]; return false;
    }
    suff = ++size; tr[cur].next[let] = size;
    tr[size].len = tr[cur].len+2;
    if (tr[size].len == 1) {
      tr[size].sufflink=2; tr[size].chain=1;
      return true;
    }
    while (true) {
      cur=tr[cur].sufflink;curlen=tr[cur].len;
     if(pos-1-curlen>=0&&s[pos-1-curlen]==s[pos]){
        tr[size].sufflink = tr[cur].next[let];
        break;
      }
    }
    tr[size].chain=1+tr[tr[size].sufflink].chain;
    return true;
  }
  void initTree() {
    memset(tr,0,sizeof tr);///CAREFUL:TESTCASES
    size = 2; suff = 2;
    tr[1].len = -1; tr[1].sufflink = 1;
    tr[2].len = 0; tr[2].sufflink = 1;
  }
}
int main() {
  int q; cin >> q;
  string operations; cin >> operations;
  PalTree::initTree();
  vector< int >subs, suffs; subs.push_back(0);
  suffs.push_back(PalTree::suff);
  for (char c : operations) {
    if (c == '-') {
      subs.pop_back(); suffs.pop_back();
      PalTree::s.pop_back();
      PalTree::suff = suffs.back();
    } else {
      PalTree::s += c;
      PalTree::addLetter(PalTree::s.size()-1);
      suffs.push_back(PalTree::suff);
      subs.push_back(subs.back()+
                PalTree::tr[PalTree::suff].chain);
    }
    cout << subs.back() << " ";
  }
	return 0;
}

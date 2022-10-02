vector<int> adj[N],cen_tree[N];
bool cent_mark[N]; int sub[N];
void dfs(int ind, int &n, int par = -1){
  n++,sub[ind] = 1;
  for(auto x : adj[ind]){
    if(x != par && !cent_mark[x])
      dfs(x, n, ind), sub[ind] += sub[x];
}}
int get_centroid(int ind,int n,int par=-1){
  for(auto x : adj[ind]){
    if(x != par && !cent_mark[x]){
      if(sub[x]>n)return get_centroid(x,n,ind);
  }} return ind;
}
int decompose(int ind){
  int n = 0; dfs(ind, n);
  int cn = get_centroid(ind, n >> 1);
  cent_mark[cn] = 1;
  for(auto x : adj[cn]){
    if(!cent_mark[x]){
      int y = decompose(x);
      cen_tree[cn].push_back(y);
  }} return cn;
}

int table[N][int(log2(N))+1], A[N];
void buildSparseforMIN(int n){
  for(int i=1;i<=n;i++) table[i][0]=A[i];
  for(int j=1; (1<<j)<=n ;j++)
    for(int i=1; (i+(1<<j)-1)<=n; i++)
      table[i][j]=min(table[i][j-1],
          table[i+(1<<(j-1))][j-1]);
}
int MIN(int i,int j){
  //int k = 32 - __builtin_clz(j-i+1) - 1;
  int k=log2(j-i+1);
  return min(table[i][k],table[j-(1<<k)+1][k]);
}

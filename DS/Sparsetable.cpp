int K[N][int(log2(N))+1], A[N];
void buildSparseforMIN(int n){
  for(int i=1;i<=n;i++) K[i][0]=A[i];
  for(int j=1; (1<<j)<=n ;j++)
    for(int i=1; (i+(1<<j)-1)<=n; i++)
      K[i][j]=min(K[i][j-1],
          K[i+(1<<(j-1))][j-1]);
}
int MIN(int i,int j){
  int k=log2(j-i+1);
  return min(K[i][k],K[j-(1<<k)+1][k]);
}

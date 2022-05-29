//dp[j][i] = min(dp[j-1][k-1] + C[k][i]) [k<=i]
//C(a,c)+C(b,d) <= C(a,d)+C(b,c) [a<=b<=c<=d]
ll dp[kmax][nmax];
void dnc(int K,int L,int R,int OptL,int OptR){
  if(L > R) return; int mid = (L + R) / 2;
  int optNow = -1; dp[K][mid] = inf;
  for(int i=OptL; i<=min(OptR, mid); i++){
    ll tmp = dp[K-1][i-1] + cost(i, mid);
    if(tmp <= dp[K][mid])
      dp[K][mid] = tmp, optNow = i;
  }
  dnc(K, L, mid - 1, OptL, optNow);
  dnc(K, mid + 1, R, optNow, OptR);
}

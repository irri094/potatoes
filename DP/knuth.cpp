//Opt[i-1][j] <= Opt[i][j] <= Opt[i][j+1]
for (int len = 2; len<=n; len++){
  for (int l=0; l+len<=n; l++){
    int r=l+len; dp[l][r] = INF;
    for(int i=opt[l][r-1]; i<=opt[l+1][r];i++){
      LL cost = dp[l][i] + dp[i][r] + C(l, r);
      if (cost < dp[l][r])
          dp[l][r] = cost, opt[l][r] = i;
    }
  }
}

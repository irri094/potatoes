/// O(segment_size*lnln(segment_size))
/// sieve generate primes upto sqrt( max high)
void segmented_sieve(ll low, ll
      high,vector<ll> &lp_of_segment) {
  int sz = high-low+1; sieve();
  for(int i=0;i<sz;i++) lp_of_segment[i]=i+low;
  for(auto p: prime) {
    if(1LL*p*p>high) break;
    for(int i = (low+p-1)/p*p-low;i<sz;i+=p) {
      if(lp_of_segment[i]==i+low)
        lp_of_segment[i] = p;
    }
  }
}

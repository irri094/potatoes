ll INF = (1ULL<<50);
struct edge { /// maintain order
	int src, dst;ll capacity;int rev; ll residue;
};
struct graph {
	int n; vector<vector<edge>> adj;
	graph(int n = 0) : n(n), adj(n) { }
	void add_edge(int src, int dst, ll capacity){
		adj[src].push_back({src, dst,
           capacity, (int)adj[dst].size()});
		adj[dst].push_back({dst, src,
           capacity, (int)adj[src].size()-1});
	}
	vector<int> level, iter;
	ll augment(int u, int t, ll cur) {
		if(u==t) return cur;
		for(int &i=iter[u]; i<adj[u].size(); ++i){
			edge &e = adj[u][i];
			if(e.residue>0 && level[u]<level[e.dst]){
				ll f=augment(e.dst, t,
               min(cur,e.residue));
				if(f>0){
          e.residue -= f;
          adj[e.dst][e.rev].residue += f;
					return f;
				}
			}
		}
		return 0;
	}
	int bfs(int s, int t) {
		level.assign(n, -1); level[s] = 0;
    queue<int> Q; Q.push(s);
		while (!Q.empty()){
			int u = Q.front(); Q.pop();
			if(u==t) break;
			for (auto &e: adj[u])
        if(e.residue>0 && level[e.dst]<0)
          Q.push(e.dst),level[e.dst]=level[u]+1;
		}
		return level[t];
	}
	ll max_flow(int s, int t){
		for(int u=0; u<n;++u) for(auto &e:adj[u])
      e.residue = e.capacity;
		ll flow = 0, itera = 0;
		while(bfs(s, t)>=0){
			iter.assign(n, 0);
			for(ll f; (f=augment(s, t, INF))>0;)
        flow += f;
		}
		return flow;
	}
	vector<edge> tree; vector<int> parent;
	void gomory_hu(){
		tree.clear(), parent.clear();
    parent.resize(n);
		for(int i=0;i<n;++i) parent[i]=0;
		for(int u = 1; u < n; ++u) {
			tree.push_back({u, parent[u],
        max_flow(u, parent[u])});
			for(int v = u+1; v < n; ++v)
        if(level[v]>=0 && parent[v]==parent[u])
          parent[v]=u;
		}
	}
};

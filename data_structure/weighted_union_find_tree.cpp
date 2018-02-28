//
// Weighted Union Find Tree
//
// verified:
//   ABC 087 D - People on a Line
//

#include <iostream>
using namespace std;

typedef int VAL;
const int MAX_UF = 110000;
const VAL SUM_UNITY = 0;

struct UnionFind {
	int par[MAX_UF];
	int rank[MAX_UF];
	VAL diff_weight[MAX_UF];

	UnionFind(int n = 0) {
		for (int i = 0; i < n; ++i) par[i] = i, rank[i] = 0, diff_weight[i] = SUM_UNITY;
	}

	void init(int n) {
		for (int i = 0; i < n; ++i) par[i] = i, rank[i] = 0, diff_weight[i] = SUM_UNITY;
	}

	int root(int x) {
		if (par[x] == x) {
			return x;
		}
		else {
			int r = root(par[x]);
			diff_weight[x] += diff_weight[par[x]];
			return par[x] = r;
		}
	}

	VAL weight(int x) {
		root(x);
		return diff_weight[x];
	}

	bool issame(int x, int y) {
		return root(x) == root(y);
	}

	bool merge(int x, int y, VAL w) {
		w += weight(x); w -= weight(y);
		x = root(x); y = root(y);
		if (x == y) return false;
		if (rank[x] < rank[y]) swap(x, y), w = -w;
		if (rank[x] == rank[y]) ++rank[x];
		par[y] = x;
		diff_weight[y] = w;
		return true;
	}

	VAL diff(int x, int y) {
		return weight(y) - weight(x);
	}

} uf;

int main() {
	int N, M;
	cin >> N >> M;
	uf.init(N);
	for (int i = 0; i < M; ++i) {
		int l, r, d;
		cin >> l >> r >> d;
		--l, --r;
		if (uf.issame(l, r)) {
			int diff = uf.diff(l, r);
			if (diff != d) {
				puts("No");
				return 0;
			}
		}
		else {
			uf.merge(l, r, d);
		}
	}
	puts("Yes");
}
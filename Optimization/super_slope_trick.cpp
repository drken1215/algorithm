//
// Super Slope Trick
//
// reference:
//   maspy: slope_super (in github)
//     https://maspypy.github.io/library/alg/monoid/add_pair.hpp
//
// verified
//   ABC 275 Ex - Monster
//     https://atcoder.jp/contests/abc275/tasks/abc275_h
//


#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")

#include <bits/stdc++.h>
using namespace std;


using ll = long long;
using u32 = uint32_t;
template<class S, class T> inline bool chmax(S &a, T b) { return (a < b ? a = b, 1 : 0); }
template<class S, class T> inline bool chmin(S &a, T b) { return (a > b ? a = b, 1 : 0); }
#define FOR1(a) for (ll _ = 0; _ < ll(a); ++_)
#define FOR2(i, a) for (ll i = 0; i < ll(a); ++i)
#define FOR3(i, a, b) for (ll i = a; i < ll(b); ++i)
#define FOR4(i, a, b, c) for (ll i = a; i < ll(b); i += (c))
#define overload4(a, b, c, d, e, ...) e
#define FOR(...) overload4(__VA_ARGS__, FOR4, FOR3, FOR2, FOR1)(__VA_ARGS__)
#define all(x) (x).begin(), (x).end()
#define len(x) ll(x.size())
#define elif else if
#define eb emplace_back
#define COUT(x) cout << #x << " = " << (x) << " (L" << __LINE__ << ")" << endl
template <class T> using vc = vector<T>;
template <class T> constexpr T infty = 0;
template <> constexpr ll infty<ll> = 2'020'000'000'000'000'000;
#define POP(a) ({ auto __x = a.back(); a.pop_back(); __x; })
template<class T> void concat(vc<T>&a, const vc<T>&b){ for(auto&x:b)a.push_back(x); }
template <typename F> ll binary_search(F check, ll ok, ll ng, bool check_ok = true) {
  if (check_ok) assert(check(ok));
  while (llabs(ok - ng) > 1) { auto x = (ng + ok) / 2; (check(x) ? ok : ng) = x; } return ok; }
template <typename E> struct Monoid_Add_Pair {
  using value_type = pair<E, E>; using X = value_type;
  static constexpr X op(const X &x, const X &y) { return {x.first + y.first, x.second + y.second}; }
  static constexpr X unit() { return {0, 0}; } static constexpr bool commute = true; };
template <class Node>
struct Node_Pool {
  struct Slot {
    union alignas(Node) {
      Slot* next;
      unsigned char storage[sizeof(Node)];
    };
  };
  using np = Node*;
  static constexpr int CHUNK_SIZE = 1 << 12;
  vc<unique_ptr<Slot[]>> chunks;
  Slot* cur = nullptr;
  int cur_used = 0;
  Slot* free_head = nullptr;
  Node_Pool() { alloc_chunk(); }
  template <class... Args>
  np create(Args&&... args) { Slot* s = new_slot(); return ::new (s) Node(forward<Args>(args)...); }
  np clone(const np x) { assert(x); Slot* s = new_slot(); return ::new (s) Node(*x); }
  void destroy(np x) { if (!x) return; x->~Node(); auto s = reinterpret_cast<Slot*>(x); s->next = free_head; free_head = s; }
  void reset() { free_head = nullptr; if (!chunks.empty()) { cur = chunks[0].get(); cur_used = 0; } }
 private:
  void alloc_chunk() { chunks.emplace_back(make_unique<Slot[]>(CHUNK_SIZE)); cur = chunks.back().get(); cur_used = 0; }
  Slot* new_slot() {
    if (free_head) { Slot* s = free_head; free_head = free_head->next; return s; }
    if (cur_used == CHUNK_SIZE) alloc_chunk();
    return &cur[cur_used++];
  }
};
template <typename Node>
struct SplayTree {
  Node_Pool<Node> pool;
  using np = Node *;
  using X = typename Node::value_type;
  using A = typename Node::operator_type;
  void free_subtree(np c) {
    if (!c) return;
    auto dfs = [&](auto &dfs, np c) -> void {
      if (c->l) dfs(dfs, c->l); if (c->r) dfs(dfs, c->r);
      c->p = c->l = c->r = nullptr; pool.destroy(c);
    };
    dfs(dfs, c);
  }
  void reset() { pool.reset(); }
  np new_root() { return nullptr; }
  np new_node(const X &x) { np n = pool.create(); Node::new_node(n, x); return n; }
  np new_node(const vc<X> &dat) {
    auto dfs = [&](auto &dfs, int l, int r) -> np {
      if (l == r) return nullptr;
      if (r == l + 1) return new_node(dat[l]);
      int m = (l + r) / 2;
      np l_root = dfs(dfs, l, m); np r_root = dfs(dfs, m + 1, r);
      np root = new_node(dat[m]); root->l = l_root, root->r = r_root;
      if (l_root) l_root->p = root; if (r_root) r_root->p = root;
      root->update(); return root;
    };
    return dfs(dfs, 0, len(dat));
  }
  u32 get_size(np root) { return (root ? root->size : 0); }
  np merge(np l_root, np r_root) {
    if (!l_root) return r_root; if (!r_root) return l_root;
    assert((!l_root->p) && (!r_root->p));
    splay_kth(r_root, 0); r_root->l = l_root; l_root->p = r_root; r_root->update(); return r_root;
  }
  np merge3(np a, np b, np c) { return merge(merge(a, b), c); }
  pair<np, np> split(np root, u32 k) {
    assert(!root || !root->p);
    if (k == 0) return {nullptr, root};
    if (k == (root->size)) return {root, nullptr};
    splay_kth(root, k - 1);
    np right = root->r; root->r = nullptr, right->p = nullptr; root->update();
    return {root, right};
  }
  tuple<np, np, np> split3(np root, u32 l, u32 r) {
    np nm, nr; tie(root, nr) = split(root, r); tie(root, nm) = split(root, l); return {root, nm, nr};
  }
  void goto_between(np &root, u32 l, u32 r) {
    if (l == 0 && r == root->size) return;
    if (l == 0) { splay_kth(root, r); root = root->l; return; }
    if (r == root->size) { splay_kth(root, l - 1); root = root->r; return; }
    splay_kth(root, r); np rp = root; root = rp->l; root->p = nullptr;
    splay_kth(root, l - 1); root->p = rp; rp->l = root; rp->update(); root = root->r;
  }
  vc<X> get_all(const np &root) {
    vc<X> res;
    auto dfs = [&](auto &dfs, np root) -> void {
      if (!root) return; root->push(); dfs(dfs, root->l); res.emplace_back(root->get()); dfs(dfs, root->r);
    };
    dfs(dfs, root); return res;
  }
  X get(np &root, u32 k) { assert(root == nullptr || !root->p); splay_kth(root, k); return root->get(); }
  void set(np &root, u32 k, const X &x) { assert(root != nullptr && !root->p); splay_kth(root, k); root->set(x); }
  X prod(np &root, u32 l, u32 r) {
    assert(root == nullptr || !root->p);
    using Mono = typename Node::Monoid_X;
    if (l == r) return Mono::unit();
    assert(0 <= l && l < r && r <= root->size);
    goto_between(root, l, r); X res = root->prod; splay(root, true); return res;
  }
  X prod(np &root) { using Mono = typename Node::Monoid_X; return (root ? root->prod : Mono::unit()); }
  void apply(np &root, u32 l, u32 r, const A &a) {
    if (l == r) return; assert(0 <= l && l < r && r <= root->size);
    goto_between(root, l, r); root->apply(a); splay(root, true);
  }
  void apply(np &root, const A &a) { if (!root) return; root->apply(a); }
  void rotate(Node *n) {
    Node *pp, *p, *c; p = n->p; pp = p->p;
    if (p->l == n) { c = n->r; n->r = p; p->l = c; } else { c = n->l; n->l = p; p->r = c; }
    if (pp && pp->l == p) pp->l = n; if (pp && pp->r == p) pp->r = n;
    n->p = pp; p->p = n; if (c) c->p = p;
  }
  void push_from_root(np c) { if (!c->p) { c->push(); return; } push_from_root(c->p); c->push(); }
  void splay(Node *me, bool push_from_root_done) {
    if (!push_from_root_done) push_from_root(me);
    me->push();
    while (me->p) {
      np p = me->p; np pp = p->p;
      if (!pp) { rotate(me); p->update(); break; }
      bool same = (p->l == me && pp->l == p) || (p->r == me && pp->r == p);
      if (same) rotate(p), rotate(me); if (!same) rotate(me), rotate(me);
      pp->update(), p->update();
    }
    me->update();
  }
  void splay_kth(np &root, u32 k) {
    assert(0 <= k && k < (root->size));
    while (1) {
      root->push();
      u32 s1 = (root->l ? root->l->size : 0);
      u32 s2 = (root->size) - (root->r ? root->r->size : 0);
      if (k < s1) root = root->l; elif (k < s2) { break; } else { k -= s2; root = root->r; }
    }
    splay(root, true);
  }
  template <typename F> pair<np, np> split_max_right(np root, F check) {
    if (!root) return {nullptr, nullptr}; assert(!root->p);
    np c = find_max_right(root, check);
    if (!c) { splay(root, true); return {nullptr, root}; }
    splay(c, true); np right = c->r; if (!right) return {c, nullptr};
    right->p = nullptr; c->r = nullptr; c->update(); return {c, right};
  }
  template <typename F> pair<np, np> split_max_right_prod(np root, F check) {
    if (!root) return {nullptr, nullptr}; assert(!root->p);
    np c = find_max_right_prod(root, check);
    if (!c) { splay(root, true); return {nullptr, root}; }
    splay(c, true); np right = c->r; if (!right) return {c, nullptr};
    right->p = nullptr; c->r = nullptr; c->update(); return {c, right};
  }
  template <typename F> np find_max_right(np root, const F &check) {
    np last_ok = nullptr, last = nullptr;
    while (root) { last = root; root->push();
      if (check(root->x)) { last_ok = root; root = root->r; } else { root = root->l; } }
    splay(last, true); return last_ok;
  }
  template <typename F> np find_max_right_prod(np root, const F &check) {
    using Mono = typename Node::Monoid_X; X prod = Mono::unit();
    np last_ok = nullptr, last = nullptr;
    while (root) { last = root; root->push();
      np tmp = root->r; root->r = nullptr; root->update();
      X lprod = Mono::op(prod, root->prod); root->r = tmp; root->update();
      if (check(lprod)) { prod = lprod; last_ok = root; root = root->r; } else { root = root->l; } }
    splay(last, true); return last_ok;
  }
};
template <typename T>
struct Node {
  using value_type = pair<T, T>; using operator_type = T; using np = Node *;
  using Monoid_X = Monoid_Add_Pair<T>;
  np p, l, r; bool rev; u32 size; pair<T, T> x; pair<T, T> prod; T add_x;
  static void new_node(np n, const pair<T, T> &x) {
    n->p = n->l = n->r = nullptr, n->rev = 0, n->size = 1;
    n->x = x, n->prod = {x.second, x.first * x.second}, n->add_x = 0;
  }
  void update() {
    size = 1; if (l) size += l->size; if (r) size += r->size;
    prod = {x.second, x.first * x.second};
    if (l) prod = Monoid_X::op(prod, l->prod);
    if (r) prod = Monoid_X::op(prod, r->prod);
  }
  void push() {
    assert(!rev); if (add_x == 0) return;
    if (l) l->x.first += add_x, l->prod.second += l->prod.first * add_x, l->add_x += add_x;
    if (r) r->x.first += add_x, r->prod.second += r->prod.first * add_x, r->add_x += add_x;
    add_x = 0;
  }
  void apply(T a) { x.first += a, prod.second += a * prod.first, add_x += a; }
  pair<T, T> get() { return x; }
  void set(const pair<T, T> &xx) { x = xx; update(); }
};
template <typename T>
struct Slope_Trick_Super {
  SplayTree<Node<T>> ST;
  using np = Node<T> *;
  struct FUNC { np root; T x0, x1, a0, y0; int size() { return (root ? root->size : 0); } };
  FUNC segment_func(T L, T R, T a, T b) { return {nullptr, L, R, a, a * L + b}; }
  FUNC from_points(vc<pair<T, T>> &point) {
    return from_points(len(point), [&](int i) -> pair<T, T> { return point[i]; }); }
  template <typename F> FUNC from_points(int N, F f) {
    vc<T> X(N), Y(N); FOR(i, N) tie(X[i], Y[i]) = f(i);
    if (N == 1) return segment_func(X[0], X[0], 0, Y[0]);
    T a0 = (Y[1] - Y[0]) / (X[1] - X[0]); T x0 = X[0], x1 = X.back();
    vc<pair<T, T>> dat; T a = a0;
    FOR(i, 1, N - 1) { T a1 = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]); dat.eb(X[i], a1 - a), a = a1; }
    return FUNC{ST.new_node(dat), x0, x1, a0, Y[0]};
  }
  pair<T, T> domain(FUNC &f) { return {f.x0, f.x1}; }
  T eval(FUNC &f, T x) {
    auto [x0, x1] = domain(f);
    if (!(x0 <= x && x <= x1)) return infty<T>;
    auto [l, r] = ST.split_max_right(f.root, [&](auto dat) -> bool { return dat.first <= x; });
    auto [a_sum, xa_sum] = ST.prod(l); f.root = ST.merge(l, r);
    return f.y0 + f.a0 * (x - x0) + a_sum * x - xa_sum;
  }
  FUNC restrict_domain(FUNC &f, T L, T R) {
    auto [x0, x1] = domain(f); chmax(L, x0), chmin(R, x1);
    if (L > R) { ST.free_subtree(f.root), f.root = nullptr; f.root = nullptr; f.x0 = infty<T>, f.x1 = -infty<T>; return f; }
    auto [l, r] = ST.split_max_right(f.root, [&](auto dat) -> bool { return dat.first < R; });
    ST.free_subtree(r);
    tie(l, r) = ST.split_max_right(l, [&](auto dat) -> bool { return dat.first <= L; });
    auto [a_sum, xa_sum] = ST.prod(l);
    T new_a0 = f.a0 + a_sum; T new_y0 = f.y0 + f.a0 * (L - x0) + a_sum * L - xa_sum;
    ST.free_subtree(l); f.root = r, f.x0 = L, f.x1 = R, f.a0 = new_a0, f.y0 = new_y0; return f;
  }

  // 定義域を「左」に広げる: [L, x1] へ拡張し、x < x0 では傾き slope の直線で延長。
  //   凸性のため slope <= f.a0 が必要（slope == f.a0 なら折れ点を作らず素直に延長）。
  //   旧左端 x0 には傾きの跳び (a0 - slope) の折れ点が入る。
  FUNC extend_left(FUNC &f, T L, T slope) {
    if (f.x0 > f.x1) return f;             // 空領域は何もしない
    if (L >= f.x0) return f;               // 拡張不要
    assert(slope <= f.a0);                 // 凸性
    T d = f.a0 - slope;                    // 旧左端 x0 での傾きの跳び
    if (d != 0) f.root = ST.merge(ST.new_node({f.x0, d}), f.root);
    f.y0 += slope * (L - f.x0);            // f(L) = f(x0) + slope*(L - x0)
    f.a0 = slope;
    f.x0 = L;
    return f;
  }
  // 定義域を「右」に広げる: [x0, R] へ拡張し、x > x1 では傾き slope の直線で延長。
  //   凸性のため slope >= (現在の右端の傾き) が必要。旧右端 x1 に折れ点が入る。
  FUNC extend_right(FUNC &f, T R, T slope) {
    if (f.x0 > f.x1) return f;
    if (R <= f.x1) return f;
    T cur_right = f.a0 + ST.prod(f.root).first;   // 現在の右端の傾き = a0 + Σ傾き変化
    assert(slope >= cur_right);            // 凸性
    T d = slope - cur_right;               // 旧右端 x1 での傾きの跳び
    if (d != 0) f.root = ST.merge(f.root, ST.new_node({f.x1, d}));
    f.x1 = R;
    return f;
  }
  // 定義域を [L,R] へ拡張（両端の傾きのまま延長。折れ点を作らず凸性を保存）。
  //   restrict_domain の逆向きに対応する版。
  FUNC extend_domain(FUNC &f, T L, T R) {
    if (f.x0 > f.x1) return f;
    if (L < f.x0) extend_left(f, L, f.a0);
    if (R > f.x1) extend_right(f, R, f.a0 + ST.prod(f.root).first);
    return f;
  }

  FUNC add(FUNC &f, FUNC &g) {
    T x0 = max(f.x0, g.x0); T x1 = min(f.x1, g.x1);
    restrict_domain(f, x0, x1), restrict_domain(g, x0, x1);
    if (x0 > x1) return f;
    T y0 = f.y0 + g.y0, a0 = f.a0 + g.a0;
    if (f.size() < g.size()) swap(f, g);
    auto tmp = ST.get_all(g.root); ST.free_subtree(g.root);
    f.x0 = x0, f.x1 = x1, f.a0 = a0, f.y0 = y0;
    if (!f.root) { f.root = ST.new_node(tmp); return f; }
    auto dfs = [&](auto &dfs, np root, int l, int r) -> void {
      if (l == r) return; root->push(); T x = root->x.first;
      int m = binary_search([&](int i) -> bool { return tmp[i].first >= x; }, r, l - 1, 0);
      if (l < m) { if (!root->l) { root->l = ST.new_node({tmp.begin() + l, tmp.begin() + m}); } else { dfs(dfs, root->l, l, m); } root->l->p = root; }
      if (m < r) { if (!root->r) { root->r = ST.new_node({tmp.begin() + m, tmp.begin() + r}); } else { dfs(dfs, root->r, m, r); } root->r->p = root; }
      root->update();
    };
    dfs(dfs, f.root, 0, len(tmp)); return f;
  }
  FUNC shift(FUNC &f, T add_x, T add_y) { ST.apply(f.root, add_x); f.x0 += add_x, f.x1 += add_x, f.y0 += add_y; return f; }
  FUNC add_const(FUNC &f, T a) { f.y0 += a; return f; }
  FUNC add_linear(FUNC &f, T a, T b) { f.y0 += a * f.x0 + b; f.a0 += a; return f; }
  pair<T, T> get_min(FUNC &f) {
    if (f.x0 > f.x1) return {infty<T>, 0};
    if (f.a0 >= 0) { return {f.y0, f.x0}; }
    auto [l, r] = ST.split_max_right_prod(f.root, [&](auto prod) -> bool { return f.a0 + prod.first < 0; });
    auto [asum, xasum] = ST.prod(l); T x = (r ? ST.get(r, 0).first : f.x1);
    f.root = ST.merge(l, r); T y = f.y0 + f.a0 * (x - f.x0) + x * asum - xasum; return {y, x};
  }
  FUNC clear_left(FUNC &f) {
    if (f.a0 >= 0) { return f; }
    auto [l, r] = ST.split_max_right_prod(f.root, [&](auto prod) -> bool { return f.a0 + prod.first < 0; });
    auto [asum, xasum] = ST.prod(l);
    if (!r) { T x = f.x1; T y = f.y0 + f.a0 * (x - f.x0) + x * asum - xasum; ST.free_subtree(l); f.root = nullptr; f.y0 = y, f.a0 = 0; return f; }
    T x = ST.get(r, 0).first; T y = f.y0 + f.a0 * (x - f.x0) + x * asum - xasum; // ★FIX: 元は ST.get(f.root,0)（split後 f.root は stale で assertion 落ち）
    T a = f.a0 + asum + ST.get(r, 0).second; ST.free_subtree(l); f.root = r; ST.set(r, 0, {x, a}); f.y0 = y; f.a0 = 0; return f;
  }
  // --- デバッグ: 折れ点座標を全部出す ---
  void dump(FUNC &f, const string& name="") {
    fprintf(stderr, "[dump %s] ", name.c_str());
    if (f.x0 > f.x1) { fprintf(stderr, "EMPTY domain\n"); return; }
    auto dat = ST.get_all(f.root);
    T x = f.x0, y = f.y0, a = f.a0;
    fprintf(stderr, "dom=[%lld,%lld] a0=%lld : ", (ll)f.x0,(ll)f.x1,(ll)f.a0);
    fprintf(stderr, "(%lld,%lld)", (ll)x,(ll)y);
    for (auto [xi, da] : dat) { y += a*(xi-x); x=xi; fprintf(stderr, " ->[slope %lld]-> (%lld,%lld)", (ll)a,(ll)x,(ll)y); a += da; }
    y += a*(f.x1-x); fprintf(stderr, " ->[slope %lld]-> (%lld,%lld)\n", (ll)a,(ll)f.x1,(ll)y);
  }
};


//------------------------------//
// Examples
//------------------------------//

// ABC 275 Ex - Monster
/*
　木を Cartesian 木（根 = 区間の最大 B）で辿り、各ノードで
　dp[v][x] := v を根とする部分木に属するモンスターについて、
             残り x 以下の攻撃を要する状態にするための消費魔力の総和の最小値

　更新時
　dp[l], dp[r] のうち x < A[v] の領域を削る
　dp[v][x] = min_{y >= x} (dp[l][y] + dp[r][y] + B[v]y) - B[v]x
　その後、x ≧ A[v] に制限された定義域を、x ≧ 0 まで広げる

　答え
　eval(dp[root], 0)
*/
void ABC_275_Ex() {
    long long N;
    cin >> N;
    vector<long long> A(N), B(N);
    for (int i = 0; i < N; i++) cin >> A[i];
    for (int i = 0; i < N; i++) cin >> B[i];
    long long MAXX = 0;
    for (auto a : A) MAXX = max(MAXX, a);

    // Cartesian 木の構築とトポロジカルソート
    vector<int> lch(N, -1), rch(N, -1), stk;
    for (int i = 0; i < N; i++) {
        int last = -1;
        while (!stk.empty() && B[stk.back()] < B[i]) { last = stk.back(); stk.pop_back(); }
        lch[i] = last;
        if (!stk.empty()) rch[stk.back()] = i;
        stk.push_back(i);
    }
    int root = stk.front();
    vector<int> order;
    vector<int> s1 = {root};
    while (!s1.empty()) {
        int v = s1.back();
        s1.pop_back();
        order.push_back(v);
        if (lch[v] != -1) s1.push_back(lch[v]);
        if (rch[v] != -1) s1.push_back(rch[v]);
    }
    reverse(order.begin(), order.end());

    Slope_Trick_Super<long long> ST;
    using F = typename decltype(ST)::FUNC;
    vector<F> dp(N);
    for (int v : order) {
        // 定義域 [0, MAXX] 上の f(x) 0
        dp[v] = ST.segment_func(0, MAXX, 0, 0);

        // 左右からのマージ
        if (lch[v] != -1) dp[v] = ST.add(dp[v], dp[lch[v]]);
        if (rch[v] != -1) dp[v] = ST.add(dp[v], dp[rch[v]]);

        // 定義域を [A[v], MAXX] に制限
        dp[v] = ST.restrict_domain(dp[v], A[v], MAXX);

        // その後の処理
        dp[v] = ST.add_linear(dp[v], B[v], 0);
        dp[v] = ST.clear_left(dp[v]);
        dp[v] = ST.add_linear(dp[v], -B[v], 0);

        // 定義域を傾き -B[v] で [0, MAXX] に拡張
        dp[v] = ST.extend_left(dp[v], 0, -B[v]);

        //COUT(v); ST.dump(dp[v]);
    }
    cout << ST.eval(dp[root], 0) << "\n";                       // 答えは eval(・, 0)
}


int main() {
    ABC_275_Ex();
}
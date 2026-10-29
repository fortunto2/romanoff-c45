// Seed histogram for the Romanoff obstruction bound, any base prime set (v2 of dp_seed:
// byte-packed keys, Booth rotations, threaded steps; same output).
//
// For base primes P0 (processed in the order given on the command line), M0 = prod q,
// T0 = lcm ord_q(2), and each a mod M0:
//   nu(a) = #{r mod T0 : a != 2^r (mod q) for every q in P0}.
// Prints "nu count" lines: #{a mod M0 : nu(a) = nu}, exact (128-bit counts).
//
// State after a prefix of primes: counts of admissible exponents per class mod L, L = lcm of
// the orders still to come. States are merged up to x -> u x + t (u a unit mod L): a prime of
// order d | L blocks the class s mod d with weight 1 for each s (weight q - d for "none"),
// and x -> u x + t permutes those classes, so the outcome distribution is equivariant and nu
// (a plain sum) invariant. The last K primes are enumerated depth-first per state.
//
// usage: dp_seed2 [--self-test] [--affine-max L] [--fuse K] [--threads N] [--split k] [--groups g q.. / q.. / ...] q1 q2 ...

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using u128 = unsigned __int128;
using Vec = std::vector<uint16_t>;
using Key = std::string;  // packed: 1 byte per entry if every entry fits, else 2
using States = std::unordered_map<Key, u128>;
struct H128 {
    size_t operator()(unsigned __int128 x) const noexcept {
        uint64_t a = uint64_t(x), b = uint64_t(x >> 64);
        a ^= b * 0x9E3779B97F4A7C15ull;
        a ^= a >> 29; a *= 0xBF58476D1CE4E5B9ull; a ^= a >> 32;
        return a;
    }
};

static int AFFINE_MAX = 420;
static int FUSE = 0;
static int THREADS = 1;
static int SPLIT = 0;  // >0: first SPLIT primes form side A
static int STREAM = 0;  // >0: last two groups streamed in this many passes
static int GROUPS_G = 0;
static bool HIER = false;  // --hier: groups after '|' are A groups (orders prime to 3)  // >0: group-product mode over Z/GROUPS_G, groups separated by '/'

static std::string u128s(u128 x) {
    if (x == 0) return "0";
    std::string s;
    while (x) { s.push_back(char('0' + int(x % 10))); x /= 10; }
    std::reverse(s.begin(), s.end());
    return s;
}

// 256-bit accumulators: the base modulus M0 can exceed 2^128 (27 primes: ~1e47), while each
// factor (an A-law or B-law weight) stays below 2^128.
struct U256 {
    u128 lo = 0, hi = 0;
    bool operator==(U256 const& o) const { return lo == o.lo && hi == o.hi; }
};
static inline void add256(U256& a, U256 const& b) {
    u128 lo = a.lo + b.lo;
    a.hi += b.hi + (lo < a.lo ? 1 : 0);
    a.lo = lo;
}
static inline U256 mul128(u128 a, u128 b) {
    const uint64_t a0 = uint64_t(a), a1 = uint64_t(a >> 64), b0 = uint64_t(b), b1 = uint64_t(b >> 64);
    const u128 p00 = (u128)a0 * b0, p01 = (u128)a0 * b1, p10 = (u128)a1 * b0, p11 = (u128)a1 * b1;
    const u128 mid = p01 + p10;
    const u128 midc = mid < p01 ? 1 : 0;
    U256 r;
    r.lo = p00 + (mid << 64);
    r.hi = p11 + (mid >> 64) + (midc << 64) + (r.lo < p00 ? 1 : 0);
    return r;
}
static inline void add_mul(U256& acc, u128 a, u128 b) { add256(acc, mul128(a, b)); }
static U256 mul_small(U256 a, uint64_t m) {
    // a * m, m < 2^64 (no overflow past 256 bits assumed)
    U256 lo = mul128(a.lo, m);
    lo.hi += a.hi * m;
    return lo;
}
static std::string u256s(U256 x) {
    if (x.hi == 0) return u128s(x.lo);
    // repeated division by 10^19 on four 64-bit limbs
    uint64_t limb[4] = {uint64_t(x.lo), uint64_t(x.lo >> 64), uint64_t(x.hi), uint64_t(x.hi >> 64)};
    const uint64_t base = 10000000000000000000ull;
    std::vector<uint64_t> parts;
    auto zero = [&] { return !limb[0] && !limb[1] && !limb[2] && !limb[3]; };
    while (!zero()) {
        u128 rem = 0;
        for (int i = 3; i >= 0; --i) {
            u128 cur = (rem << 64) | limb[i];
            limb[i] = uint64_t(cur / base);
            rem = cur % base;
        }
        parts.push_back(uint64_t(rem));
    }
    std::string s = std::to_string(parts.back());
    for (int i = int(parts.size()) - 2; i >= 0; --i) {
        std::string p = std::to_string(parts[i]);
        s += std::string(19 - p.size(), '0') + p;
    }
    return s;
}

static int ord2(int q) {
    int x = 2 % q, d = 1;
    while (x != 1) { x = int((int64_t)x * 2 % q); ++d; }
    return d;
}

static Key pack(Vec const& v, int width) {
    Key k(v.size() * width, '\0');
    for (size_t i = 0; i < v.size(); ++i) {
        if (width == 1) k[i] = char(v[i]);
        else { k[2 * i] = char(v[i] >> 8); k[2 * i + 1] = char(v[i] & 255); }
    }
    return k;
}

static void unpack(Key const& k, int width, Vec& v) {
    v.resize(k.size() / width);
    for (size_t i = 0; i < v.size(); ++i)
        v[i] = width == 1 ? uint8_t(k[i]) : uint16_t(uint8_t(k[2 * i]) << 8 | uint8_t(k[2 * i + 1]));
}

// Booth: start of the lexicographically least rotation, O(L).
static int least_rotation(uint16_t const* s, int n, std::vector<int>& f) {
    f.assign(2 * n, -1);
    int k = 0;
    for (int j = 1; j < 2 * n; ++j) {
        uint16_t sj = s[j % n];
        int i = f[j - k - 1];
        while (i != -1 && sj != s[(k + i + 1) % n]) {
            if (sj < s[(k + i + 1) % n]) k = j - i - 1;
            i = f[i];
        }
        if (i == -1 && sj != s[(k + i + 1) % n]) {
            if (sj < s[(k + i + 1) % n]) k = j;
            f[j - k] = -1;
        } else {
            f[j - k] = i + 1;
        }
    }
    return k;
}

struct Units {
    std::map<int, std::vector<std::vector<int>>> perms;
    std::vector<std::vector<int>> const& get(int L) {
        auto it = perms.find(L);
        if (it != perms.end()) return it->second;
        std::vector<std::vector<int>> ps;
        for (int u = 1; u < L; ++u) {
            if (std::gcd(u, L) != 1) continue;
            std::vector<int> p(L);
            for (int x = 0; x < L; ++x) p[x] = int((int64_t)u * x % L);
            ps.push_back(std::move(p));
        }
        return perms[L] = std::move(ps);
    }
};

struct Canon {
    std::vector<std::vector<int>> const* ps = nullptr;
    Vec w, r, best;
    std::vector<int> f;
    int L = 0;
    explicit Canon(int L_, Units& U) : L(L_) {
        if (L > 2 && L <= AFFINE_MAX) ps = &U.get(L);
        w.resize(L); r.resize(L);
    }
    void rot(Vec const& a, Vec& out) {
        int k = least_rotation(a.data(), L, f);
        for (int i = 0; i < L; ++i) out[i] = a[(i + k) % L];
    }
    // canonical form of v, written into best
    Vec const& operator()(Vec const& v) {
        best.resize(L);
        if (L <= 1) { best = v; return best; }
        if (!ps) { rot(v, best); return best; }
        bool first = true;
        for (auto const& p : *ps) {
            for (int x = 0; x < L; ++x) w[p[x]] = v[x];
            rot(w, r);
            if (first || r < best) { best = r; first = false; }
        }
        return best;
    }
};

struct Hist {
    std::map<int, u128> h;                // g == 1: nu -> count
    std::map<int, U256> h256;             // hier mode: exact counts when M0 >= 2^128
    U256 M256;
    std::map<Vec, u128> fin;              // g > 1: canonical final vector over Z/g -> count
    u128 M = 1;
    int T = 1;
    int g = 1;
};

static Hist run(std::vector<int> const& primes, bool trace, int g = 1) {
    const size_t n = primes.size();
    std::vector<int> d(n);
    Hist out;
    for (size_t i = 0; i < n; ++i) {
        d[i] = ord2(primes[i]);
        out.M *= (u128)primes[i];
        out.T = std::lcm(out.T, d[i]);
    }
    if (out.T > 65535) throw std::runtime_error("T0 too large for uint16 counts");
    std::vector<int> L(n + 1, g);
    for (size_t i = n; i-- > 0;) L[i] = std::lcm(L[i + 1], d[i]);
    out.T = L[0];
    out.g = g;
    auto width = [&](size_t i) { return out.T / L[i] <= 255 ? 1 : 2; };

    Units U;
    States st;
    {
        Canon c(L[0], U);
        st[pack(c(Vec(L[0], uint16_t(out.T / L[0]))), width(0))] = 1;
    }
    if (trace) std::cerr << "M0 = " << u128s(out.M) << "\nT0 = " << out.T << "\n";
    const size_t K = std::min<size_t>(FUSE, n);
    for (size_t i = 0; i + K < n; ++i) {
        const int q = primes[i], dq = d[i], Lc = L[i], Ln = L[i + 1];
        const int wc = width(i), wn = width(i + 1);
        U.get(Ln <= AFFINE_MAX && Ln > 2 ? Ln : 3);  // build tables before threads start
        U.get(g <= AFFINE_MAX && g > 2 ? g : 3);
        std::vector<States> part(THREADS);
        auto worker = [&](int tid) {
            Canon c(Ln, U);
            Vec v, o(Ln);
            States& nx = part[tid];
            size_t idx = 0;
            for (auto const& [k, w] : st) {
                if (idx++ % THREADS != (size_t)tid) continue;
                unpack(k, wc, v);
                std::fill(o.begin(), o.end(), 0);
                for (int x = 0; x < Lc; ++x) o[x % Ln] += v[x];
                nx[pack(c(o), wn)] += w * (u128)(q - dq);
                for (int s = 0; s < dq; ++s) {
                    std::fill(o.begin(), o.end(), 0);
                    for (int x = 0; x < Lc; ++x)
                        if (x % dq != s) o[x % Ln] += v[x];
                    nx[pack(c(o), wn)] += w;
                }
            }
        };
        if (THREADS == 1) worker(0);
        else {
            std::vector<std::thread> ts;
            for (int t = 0; t < THREADS; ++t) ts.emplace_back(worker, t);
            for (auto& t : ts) t.join();
        }
        st.clear();
        st.swap(part[0]);
        for (int t = 1; t < THREADS; ++t) {
            for (auto& [k, w] : part[t]) st[k] += w;
            States().swap(part[t]);
        }
        if (trace)
            std::cerr << "step " << i + 1 << ": q=" << q << " ord=" << dq << " L=" << Lc
                      << "->" << Ln << " states=" << st.size() << std::endl;
    }

    u128 tot = 0;
    const size_t i0 = n - K;
    std::vector<std::vector<u128>> hs(THREADS, std::vector<u128>(out.T + 1, 0));
    std::vector<std::map<Vec, u128>> fins(THREADS);
    // small g: leaves go to a raw hash map on a packed key, canonicalised once at the end
    const bool packed = g > 1 && g <= 16 && out.T / g <= 255;
    std::vector<std::unordered_map<u128, u128, H128>> raw(THREADS);
    auto tail = [&](int tid) {
        Canon cg(g, U);
        std::vector<Vec> buf(n + 1);
        for (size_t i = i0; i <= n; ++i) buf[i].assign(L[i], 0);
        auto& h = hs[tid];
        std::function<void(size_t, Vec const&, u128)> dfs = [&](size_t i, Vec const& v, u128 w) {
            if (i == n) {
                if (g == 1) h[v[0]] += w;
                else if (packed) {
                    u128 key = 0;
                    for (int c = 0; c < g; ++c) key |= u128(v[c]) << (8 * c);
                    raw[tid][key] += w;
                } else fins[tid][cg(v)] += w;
                return;
            }
            const int q = primes[i], dq = d[i], Lc = L[i], Ln = L[i + 1];
            if (g == 1 && i + 1 == n) {  // Lc == dq: nu = sum - v[s]
                int S = 0;
                for (int x = 0; x < Lc; ++x) S += v[x];
                h[S] += w * (u128)(q - dq);
                for (int s = 0; s < dq; ++s) h[S - v[s]] += w;
                return;
            }
            Vec& o = buf[i + 1];
            std::fill(o.begin(), o.end(), 0);
            for (int x = 0; x < Lc; ++x) o[x % Ln] += v[x];
            dfs(i + 1, o, w * (u128)(q - dq));
            for (int s = 0; s < dq; ++s) {
                std::fill(o.begin(), o.end(), 0);
                for (int x = 0; x < Lc; ++x)
                    if (x % dq != s) o[x % Ln] += v[x];
                dfs(i + 1, o, w);
            }
        };
        Vec v;
        size_t idx = 0;
        for (auto const& [k, w] : st) {
            if (idx++ % THREADS != (size_t)tid) continue;
            unpack(k, width(i0), v);
            dfs(i0, v, w);
        }
    };
    if (THREADS == 1) tail(0);
    else {
        std::vector<std::thread> ts;
        for (int t = 0; t < THREADS; ++t) ts.emplace_back(tail, t);
        for (auto& t : ts) t.join();
    }
    for (int nu = 0; nu <= out.T; ++nu) {
        u128 c = 0;
        for (auto const& h : hs) c += h[nu];
        if (c) { out.h[nu] = c; tot += c; }
    }
    for (auto& f : fins)
        for (auto const& [v, w] : f) { out.fin[v] += w; tot += w; }
    if (packed) {
        Canon cg(g, U);
        Vec v(g);
        for (auto& r : raw) {
            for (auto const& [key, w] : r) {
                for (int c = 0; c < g; ++c) v[c] = uint16_t((key >> (8 * c)) & 255);
                out.fin[cg(v)] += w;
                tot += w;
            }
            std::unordered_map<u128, u128, H128>().swap(r);
        }
    }
    if (trace && K) std::cerr << "tail: " << K << " primes enumerated from " << st.size() << " states" << std::endl;
    if (tot != out.M) throw std::runtime_error("histogram total != M0");
    return out;
}

// Split the base into A (first k primes) and B. With L_A, L_B their periods, g = gcd and
// T0 = lcm, r mod T0 <-> (r mod L_A, r mod L_B) agreeing mod g, so
//   nu = sum_{c mod g} alpha_A(c) alpha_B(c),  alpha_X(c) = #{x mod L_X : x = c (g), admissible}.
// Each side's distribution of alpha is invariant under the group it was canonicalised by
// (affine maps of Z/g lift to Z/L_X), so pairing orbit representatives of A (total orbit
// weight) with every element of B's orbits (per-element weight) gives the exact histogram.
static Hist split_run(std::vector<int> const& A, std::vector<int> const& B, bool trace) {
    int LA = 1, LB = 1;
    for (int q : A) LA = std::lcm(LA, ord2(q));
    for (int q : B) LB = std::lcm(LB, ord2(q));
    const int g = std::gcd(LA, LB);
    if (trace) std::cerr << "split: L_A=" << LA << " L_B=" << LB << " g=" << g << "\n";
    Hist ha = run(A, trace, g), hb = run(B, trace, g);
    Hist out;
    out.M = ha.M * hb.M;
    out.T = std::lcm(LA, LB);
    if (out.T > 65535) throw std::runtime_error("T0 too large");
    if (g == 1) {  // independent: nu = nu_A * nu_B
        for (auto const& [a, wa] : ha.h) for (auto const& [b, wb] : hb.h) out.h[a * b] += wa * wb;
        return out;
    }
    // group used by Canon(g): translations, times units when 2 < g <= AFFINE_MAX
    std::vector<int> units{1};
    if (g > 2 && g <= AFFINE_MAX) {
        units.clear();
        for (int u = 1; u < g; ++u) if (std::gcd(u, g) == 1) units.push_back(u);
    }
    auto orbit_size = [&](Vec const& b) {
        std::vector<Vec> orbit;
        for (int u : units)
            for (int t = 0; t < g; ++t) {
                Vec im(g);
                for (int x = 0; x < g; ++x) im[(u * x + t) % g] = b[x];
                orbit.push_back(std::move(im));
            }
        std::sort(orbit.begin(), orbit.end());
        return size_t(std::unique(orbit.begin(), orbit.end()) - orbit.begin());
    };
    {
        size_t ea = 0, eb = 0;
        for (auto const& [a, w] : ha.fin) ea += orbit_size(a);
        for (auto const& [b, w] : hb.fin) eb += orbit_size(b);
        if (trace)
            std::cerr << "A: " << ha.fin.size() << " orbits, " << ea << " elements; B: " << hb.fin.size()
                      << " orbits, " << eb << " elements\n";
        if (ha.fin.size() * eb > hb.fin.size() * ea) std::swap(ha, hb);  // reps from the side that saves work
    }
    std::vector<std::pair<Vec, u128>> belems;
    for (auto const& [b, wb] : hb.fin) {
        std::vector<Vec> orbit;
        for (int u : units)
            for (int t = 0; t < g; ++t) {
                Vec im(g);
                for (int x = 0; x < g; ++x) im[(u * x + t) % g] = b[x];
                orbit.push_back(std::move(im));
            }
        std::sort(orbit.begin(), orbit.end());
        orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
        if (wb % (u128)orbit.size()) throw std::runtime_error("orbit weight not divisible");
        for (auto& e : orbit) belems.emplace_back(std::move(e), wb / (u128)orbit.size());
    }
    if (trace)
        std::cerr << "pairing " << ha.fin.size() << " A-orbits x " << belems.size() << " B-elements\n";
    std::vector<std::vector<u128>> hs(THREADS, std::vector<u128>(out.T + 1, 0));
    std::vector<std::pair<Vec, u128>> aorb(ha.fin.begin(), ha.fin.end());
    auto worker = [&](int tid) {
        auto& h = hs[tid];
        for (size_t i = tid; i < aorb.size(); i += THREADS) {
            Vec const& a = aorb[i].first;
            const u128 wa = aorb[i].second;
            for (auto const& [b, wb] : belems) {
                int nu = 0;
                for (int c = 0; c < g; ++c) nu += a[c] * b[c];
                h[nu] += wa * wb;
            }
        }
    };
    std::vector<std::thread> ts;
    for (int t = 0; t < THREADS; ++t) ts.emplace_back(worker, t);
    for (auto& t : ts) t.join();
    u128 tot = 0;
    for (int nu = 0; nu <= out.T; ++nu) {
        u128 c = 0;
        for (auto const& h : hs) c += h[nu];
        if (c) { out.h[nu] = c; tot += c; }
    }
    if (tot != out.M) throw std::runtime_error("split total != M0");
    return out;
}

// Group product over an interface Z/g. Groups G_1..G_K with periods L_G (lcm of their orders
// and g) and private parts P_G = L_G/g. If T0 = g * prod P_G, then r mod T0 <-> (r mod g,
// private coordinates) and every group sees only its own private coordinate, so
//   nu = sum_{c mod g} prod_G alpha_G(c),   alpha_G(c) = #{x mod L_G : x = c (g), admissible}
// with the alpha_G independent. Partial products are merged up to the affine maps of Z/g
// (the joint law is invariant under simultaneous maps); the last group is paired with them.
using GV = std::array<uint16_t, 16>;
struct GVHash {
    size_t operator()(GV const& v) const noexcept {
        uint64_t h = 1469598103934665603ull;
        for (auto x : v) { h ^= x; h *= 1099511628211ull; }
        return h;
    }
};
using GMap = std::unordered_map<GV, u128, GVHash>;
static size_t RAW_LIMIT = size_t(1) << 23;
static size_t SIGMA_LIMIT = size_t(1) << 24;

static std::vector<std::vector<int>> group_perms(int g) {
    std::vector<int> units{1};
    if (g > 2 && g <= AFFINE_MAX) {
        units.clear();
        for (int u = 1; u < g; ++u) if (std::gcd(u, g) == 1) units.push_back(u);
    }
    std::vector<std::vector<int>> ps;
    for (int u : units)
        for (int t = 0; t < g; ++t) {
            std::vector<int> p(g);
            for (int x = 0; x < g; ++x) p[x] = (u * x + t) % g;
            ps.push_back(std::move(p));
        }
    return ps;
}

static GV gcanon(GV const& v, int g, std::vector<std::vector<int>> const& ps) {
    GV best{}, w{};
    bool first = true;
    for (auto const& p : ps) {
        for (int x = 0; x < g; ++x) w[p[x]] = v[x];
        if (first || std::lexicographical_compare(w.begin(), w.begin() + g, best.begin(), best.begin() + g)) {
            best = w;
            first = false;
        }
    }
    return best;
}

static Hist group_run(std::vector<std::vector<int>> const& groups, int g, bool trace) {
    if (g > 16) throw std::runtime_error("interface g > 16");
    const auto ps = group_perms(g);
    std::vector<std::vector<std::pair<GV, u128>>> elems;
    Hist out;
    out.T = g;
    long long privprod = 1;
    for (auto const& G : groups) {
        Hist hg = run(G, false, g);
        out.M *= hg.M;
        out.T = std::lcm(out.T, hg.T);
        privprod *= hg.T / g;
        std::vector<std::pair<GV, u128>> el;
        u128 tot = 0;
        for (auto const& [b, wb] : hg.fin) {
            std::vector<GV> orbit;
            GV bb{};
            for (int x = 0; x < g; ++x) bb[x] = b[x];
            for (auto const& p : ps) {
                GV im{};
                for (int x = 0; x < g; ++x) im[p[x]] = bb[x];
                orbit.push_back(im);
            }
            std::sort(orbit.begin(), orbit.end());
            orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
            if (wb % (u128)orbit.size()) throw std::runtime_error("orbit weight not divisible");
            for (auto& e : orbit) el.emplace_back(e, wb / (u128)orbit.size());
            tot += wb;
        }
        if (tot != hg.M) throw std::runtime_error("group weight mismatch");
        if (trace) std::cerr << "group L=" << hg.T << " primes=" << G.size() << ": " << hg.fin.size()
                             << " orbits, " << el.size() << " elements" << std::endl;
        elems.push_back(std::move(el));
    }
    if ((long long)out.T != g * privprod)
        throw std::runtime_error("private parts not independent: T0 != g * prod(L_G/g)");
    if (out.T / g > 65535) throw std::runtime_error("T0/g too large for uint16 products");
    const size_t K = groups.size();
    GMap st;
    for (auto const& [v, w] : elems[0]) st[gcanon(v, g, ps)] += w;
    const bool stream = STREAM > 0 && K >= 3;
    for (size_t k = 1; k + 1 + (stream ? 1 : 0) < K; ++k) {
        std::vector<GMap> outs(THREADS);
        std::vector<std::pair<GV, u128>> sv(st.begin(), st.end());
        GMap().swap(st);
        auto worker = [&](int tid) {
            GMap raw;
            GMap& o = outs[tid];
            auto flush = [&] {
                for (auto const& [v, w] : raw) o[gcanon(v, g, ps)] += w;
                GMap().swap(raw);
            };
            for (size_t i = tid; i < sv.size(); i += THREADS) {
                auto const& [a, wa] = sv[i];
                for (auto const& [b, wb] : elems[k]) {
                    GV p{};
                    for (int c = 0; c < g; ++c) p[c] = uint16_t(a[c] * b[c]);
                    raw[p] += wa * wb;
                }
                if (raw.size() > RAW_LIMIT) flush();
            }
            flush();
        };
        std::vector<std::thread> ts;
        for (int t = 0; t < THREADS; ++t) ts.emplace_back(worker, t);
        for (auto& t : ts) t.join();
        st.swap(outs[0]);
        for (int t = 1; t < THREADS; ++t) {
            for (auto const& [v, w] : outs[t]) st[v] += w;
            GMap().swap(outs[t]);
        }
        if (trace) std::cerr << "product with group " << k + 1 << ": states=" << st.size() << std::endl;
    }
    std::vector<std::vector<u128>> hs(THREADS, std::vector<u128>(out.T + 1, 0));
    std::vector<std::pair<GV, u128>> sv(st.begin(), st.end());
    auto const& last = elems[K - 1];
    auto pairw = [&](int tid) {
        auto& h = hs[tid];
        for (size_t i = tid; i < sv.size(); i += THREADS) {
            auto const& [a, wa] = sv[i];
            for (auto const& [b, wb] : last) {
                int nu = 0;
                for (int c = 0; c < g; ++c) nu += a[c] * b[c];
                h[nu] += wa * wb;
            }
        }
    };
    if (stream) {
        // Last two groups without holding the merged states: products with group K-1 are
        // split into STREAM passes x THREADS shards by an orbit-invariant hash (sum of a
        // random table over the entries), so each shard canonicalises and merges only its
        // own orbits, then pairs them with group K. Pairing is linear, so shards add up.
        auto const& pen = elems[K - 2];
        const uint64_t shards = uint64_t(STREAM) * THREADS;
        // structure-of-arrays copies for the vector loops
        const size_t NP = pen.size(), NL = last.size();
        std::vector<std::vector<uint16_t>> PB(g, std::vector<uint16_t>(NP)), LB(g, std::vector<uint16_t>(NL));
        for (size_t j = 0; j < NP; ++j) for (int c = 0; c < g; ++c) PB[c][j] = pen[j].first[c];
        for (size_t j = 0; j < NL; ++j) for (int c = 0; c < g; ++c) LB[c][j] = last[j].first[c];
        // last-group weights as u64 when their total fits (local histogram in u64)
        u128 wl_tot = 0;
        for (auto const& e : last) wl_tot += e.second;
        const bool fast = wl_tot < (u128(1) << 63);
        std::vector<uint64_t> lw(NL);
        for (size_t j = 0; j < NL; ++j) lw[j] = uint64_t(last[j].second);
        size_t total_states = 0;
        std::mutex mu;
        for (int pass = 0; pass < STREAM; ++pass) {
            auto shard = [&](int tid) {
                const uint64_t me = uint64_t(pass) * THREADS + tid;
                GMap m;
                std::vector<uint64_t> hv(NP);
                for (auto const& [a, wa] : sv) {
                    // orbit-invariant hash: symmetric in the entries (sum of a mixed function)
                    std::fill(hv.begin(), hv.end(), 0);
                    for (int c = 0; c < g; ++c) {
                        const uint64_t ac = a[c];
                        uint16_t const* pb = PB[c].data();
                        for (size_t j = 0; j < NP; ++j) {
                            uint64_t x = ac * pb[j] + 1;
                            hv[j] += x * x * 0x9E3779B97F4A7C15ull + x * 0xC2B2AE3D27D4EB4Full;
                        }
                    }
                    for (size_t j = 0; j < NP; ++j) {
                        if ((hv[j] >> 21) % shards != me) continue;
                        GV p{};
                        for (int c = 0; c < g; ++c) p[c] = uint16_t(a[c] * pen[j].first[c]);
                        m[gcanon(p, g, ps)] += wa * pen[j].second;
                    }
                }
                auto& hh = hs[tid];
                std::vector<int> nu(NL);
                std::vector<uint64_t> loc(out.T + 1, 0);
                std::vector<int> touched;
                touched.reserve(NL);
                for (auto const& [v, w] : m) {
                    std::fill(nu.begin(), nu.end(), 0);
                    for (int c = 0; c < g; ++c) {
                        const int vc = v[c];
                        uint16_t const* lb = LB[c].data();
                        for (size_t j = 0; j < NL; ++j) nu[j] += vc * lb[j];
                    }
                    if (fast) {
                        for (size_t j = 0; j < NL; ++j) {
                            if (!loc[nu[j]]) touched.push_back(nu[j]);
                            loc[nu[j]] += lw[j];
                        }
                        for (int t : touched) { hh[t] += w * (u128)loc[t]; loc[t] = 0; }
                        touched.clear();
                    } else {
                        for (size_t j = 0; j < NL; ++j) hh[nu[j]] += w * last[j].second;
                    }
                }
                std::lock_guard<std::mutex> lk(mu);
                total_states += m.size();
            };
            std::vector<std::thread> ts;
            for (int t = 0; t < THREADS; ++t) ts.emplace_back(shard, t);
            for (auto& t : ts) t.join();
            if (trace) std::cerr << "stream pass " << pass + 1 << "/" << STREAM << ": states so far " << total_states << std::endl;
        }
    } else if (K == 1) {
        for (auto const& [v, w] : st) {
            int nu = 0;
            for (int c = 0; c < g; ++c) nu += v[c];
            hs[0][nu] += w;
        }
    } else {
        std::vector<std::thread> ts;
        for (int t = 0; t < THREADS; ++t) ts.emplace_back(pairw, t);
        for (auto& t : ts) t.join();
    }
    u128 tot = 0;
    for (int nu = 0; nu <= out.T; ++nu) {
        u128 c = 0;
        for (auto const& h : hs) c += h[nu];
        if (c) { out.h[nu] = c; tot += c; }
    }
    if (tot != out.M) throw std::runtime_error("group product total != M0");
    if (trace) std::cerr << "pairing: " << sv.size() << " states x " << last.size() << " elements" << std::endl;
    return out;
}

// Two-level interface. B groups meet in r mod 12 (as in group_run); A groups have orders
// prime to 3, so their alpha depends on c mod 4 only. With sigma_s = sum_{c = s (4)} pi_B(c),
//   nu = sum_{s mod 4} sigma_s * beta_A(s),   beta_A = prod over A groups of alpha (on Z/4).
// The law of sigma is invariant under the affine maps of Z/4 (they lift to Z/12), so sigma
// orbit representatives (total orbit weight) are paired with every element of the A law.
static std::vector<std::pair<GV, u128>> expand_orbits(Hist const& hg, int g,
                                                      std::vector<std::vector<int>> const& ps) {
    std::vector<std::pair<GV, u128>> el;
    u128 tot = 0;
    for (auto const& [b, wb] : hg.fin) {
        std::vector<GV> orbit;
        GV bb{};
        for (int x = 0; x < g; ++x) bb[x] = b[x];
        for (auto const& p : ps) {
            GV im{};
            for (int x = 0; x < g; ++x) im[p[x]] = bb[x];
            orbit.push_back(im);
        }
        std::sort(orbit.begin(), orbit.end());
        orbit.erase(std::unique(orbit.begin(), orbit.end()), orbit.end());
        if (wb % (u128)orbit.size()) throw std::runtime_error("orbit weight not divisible");
        for (auto& e : orbit) el.emplace_back(e, wb / (u128)orbit.size());
        tot += wb;
    }
    if (tot != hg.M) throw std::runtime_error("group weight mismatch");
    return el;
}

static Hist hier_run(std::vector<std::vector<int>> const& Bg, std::vector<std::vector<int>> const& Ag,
                     bool trace) {
    const int g = 12, h = 4;
    const auto ps = group_perms(g), ps4 = group_perms(h);
    Hist out;
    out.T = g;
    long long privprod = 1;
    std::vector<std::vector<std::pair<GV, u128>>> be;
    for (auto const& G : Bg) {
        Hist hg = run(G, false, g);
        out.M *= hg.M;
        out.T = std::lcm(out.T, hg.T);
        privprod *= hg.T / g;
        be.push_back(expand_orbits(hg, g, ps));
        if (trace) std::cerr << "B group L=" << hg.T << ": " << be.back().size() << " elements" << std::endl;
    }
    // A law: product of A groups on Z/4, full elements (raw merge)
    std::unordered_map<uint64_t, u128> alaw{{0x0001000100010001ull, 1}};
    long long aprod = 1;
    for (auto const& G : Ag) {
        Hist hg = run(G, false, h);
        if (hg.T % 3 == 0) throw std::runtime_error("A group order divisible by 3");
        out.M *= hg.M;
        out.T = std::lcm(out.T, hg.T);
        privprod *= hg.T / h;
        aprod *= hg.T / h;
        if (aprod > 65535) throw std::runtime_error("A entries overflow uint16");
        auto el = expand_orbits(hg, h, ps4);
        std::unordered_map<uint64_t, u128> nx;
        for (auto const& [k, w] : alaw)
            for (auto const& [b, wb] : el) {
                uint64_t key = 0;
                for (int s = 0; s < h; ++s) key |= uint64_t(((k >> (16 * s)) & 0xFFFF) * b[s]) << (16 * s);
                nx[key] += w * wb;
            }
        alaw.swap(nx);
        if (trace) std::cerr << "A group L=" << hg.T << ": " << el.size() << " elements, A law " << alaw.size() << std::endl;
    }
    if ((long long)out.T != g * privprod) throw std::runtime_error("private parts not independent");
    if (out.T / g > 65535 || out.T > 2000000000) throw std::runtime_error("T0 too large");
    // B: in-memory product of all but the last B group
    const size_t K = Bg.size();
    GMap st;
    for (auto const& [v, w] : be[0]) st[gcanon(v, g, ps)] += w;
    for (size_t k = 1; k + 1 < K; ++k) {
        GMap nx;
        for (auto const& [a, wa] : st)
            for (auto const& [b, wb] : be[k]) {
                GV p{};
                for (int c = 0; c < g; ++c) p[c] = uint16_t(a[c] * b[c]);
                nx[gcanon(p, g, ps)] += wa * wb;
            }
        st.swap(nx);
        if (trace) std::cerr << "B product " << k + 1 << ": states=" << st.size() << std::endl;
    }
    // B last group: fold the products to sigma on Z/4, canonical under Aff(Z/4)
    std::vector<std::pair<GV, u128>> sv(st.begin(), st.end());
    GMap().swap(st);
    auto const& last = be[K - 1];
    const size_t NL = last.size();
    std::vector<std::vector<uint16_t>> LB(g, std::vector<uint16_t>(NL));
    for (size_t j = 0; j < NL; ++j) for (int c = 0; c < g; ++c) LB[c][j] = last[j].first[c];
    auto canon4 = [&](uint64_t key) {
        uint64_t best = ~0ull;
        for (auto const& p : ps4) {
            uint64_t k2 = 0;
            for (int s = 0; s < h; ++s) k2 |= ((key >> (16 * s)) & 0xFFFF) << (16 * p[s]);
            // lexicographic order on (s0, s1, s2, s3): compare as big-endian
            uint64_t be2 = 0;
            for (int s = 0; s < h; ++s) be2 = (be2 << 16) | ((k2 >> (16 * s)) & 0xFFFF);
            best = std::min(best, be2);
        }
        uint64_t k3 = 0;  // back to little-endian slots
        for (int s = 0; s < h; ++s) k3 |= ((best >> (16 * (h - 1 - s))) & 0xFFFF) << (16 * s);
        return k3;
    };
    // A law as arrays for the pairing
    std::vector<std::pair<uint64_t, u128>> avec(alaw.begin(), alaw.end());
    std::unordered_map<uint64_t, u128>().swap(alaw);
    // order by entry sum so that, for a fixed sigma, nu moves slowly along the array
    auto esum = [&](uint64_t k) { uint64_t s = 0; for (int t = 0; t < h; ++t) s += (k >> (16 * t)) & 0xFFFF; return s; };
    std::sort(avec.begin(), avec.end(), [&](auto const& x, auto const& y) {
        uint64_t a = esum(x.first), b = esum(y.first);
        return a != b ? a < b : x.first < y.first;
    });
    const size_t NA = avec.size();
    std::vector<std::vector<uint32_t>> AB(h, std::vector<uint32_t>(NA));
    for (size_t j = 0; j < NA; ++j) for (int s = 0; s < h; ++s) AB[s][j] = uint32_t((avec[j].first >> (16 * s)) & 0xFFFF);
    std::vector<std::vector<U256>> hs(THREADS, std::vector<U256>(out.T + 1));
    std::vector<size_t> paired(THREADS, 0), flushes(THREADS, 0);
    // sigma orbits are paired with the A law and dropped whenever the map gets large: pairing
    // is linear, so an orbit paired in two flushes with partial weights is still exact
    auto pair_out = [&](int tid, std::unordered_map<uint64_t, u128>& o) {
        auto& hh = hs[tid];
        std::vector<uint64_t> nu(NA);
        std::vector<u128> loc(out.T + 1, 0);  // A weights for one sigma; sum <= M_A < 2^128
        std::vector<uint32_t> touched;
        touched.reserve(NA);
        for (auto const& [key, w] : o) {
            std::fill(nu.begin(), nu.end(), 0);
            for (int s = 0; s < h; ++s) {
                const uint64_t vs = (key >> (16 * s)) & 0xFFFF;
                uint32_t const* ab = AB[s].data();
                for (size_t j = 0; j < NA; ++j) nu[j] += vs * ab[j];
            }
            for (size_t j = 0; j < NA; ++j) {
                if (!loc[nu[j]]) touched.push_back(uint32_t(nu[j]));
                loc[nu[j]] += avec[j].second;
            }
            for (uint32_t t : touched) { add_mul(hh[t], w, loc[t]); loc[t] = 0; }
            touched.clear();
        }
        paired[tid] += o.size();
        flushes[tid] += 1;
        std::unordered_map<uint64_t, u128>().swap(o);
    };
    auto fold = [&](int tid) {
        std::unordered_map<uint64_t, u128> raw, o;
        std::vector<uint32_t> part(h * NL);
        for (size_t i = tid; i < sv.size(); i += THREADS) {
            auto const& [a, wa] = sv[i];
            std::fill(part.begin(), part.end(), 0);
            for (int c = 0; c < g; ++c) {
                const uint32_t ac = a[c];
                uint16_t const* lb = LB[c].data();
                uint32_t* pp = part.data() + size_t(c % h) * NL;
                for (size_t j = 0; j < NL; ++j) pp[j] += ac * lb[j];
            }
            for (size_t j = 0; j < NL; ++j) {
                uint64_t key = 0;
                for (int s = 0; s < h; ++s) key |= uint64_t(part[s * NL + j]) << (16 * s);
                raw[key] += wa * last[j].second;
            }
            if (raw.size() > RAW_LIMIT) {
                for (auto const& [k, w] : raw) o[canon4(k)] += w;
                std::unordered_map<uint64_t, u128>().swap(raw);
                if (o.size() > SIGMA_LIMIT) pair_out(tid, o);
            }
        }
        for (auto const& [k, w] : raw) o[canon4(k)] += w;
        std::unordered_map<uint64_t, u128>().swap(raw);
        pair_out(tid, o);
    };
    {
        std::vector<std::thread> ts;
        for (int t = 0; t < THREADS; ++t) ts.emplace_back(fold, t);
        for (auto& t : ts) t.join();
    }
    if (trace) {
        size_t P = 0, F = 0;
        for (int t = 0; t < THREADS; ++t) { P += paired[t]; F += flushes[t]; }
        std::cerr << "sigma orbits paired: " << P << " in " << F << " flushes; A law " << NA << std::endl;
    }
    U256 tot;
    for (size_t nu = 0; nu <= size_t(out.T); ++nu) {
        U256 c;
        for (auto const& hh : hs) add256(c, hh[nu]);
        if (c.lo || c.hi) {
            out.h256[int(nu)] = c;
            out.h[int(nu)] = c.lo;  // exact only while M0 < 2^128 (self-tests)
            add256(tot, c);
        }
    }
    U256 M0;
    M0.lo = 1;
    for (auto const* G : {&Bg, &Ag})
        for (auto const& gr : *G)
            for (int q : gr) M0 = mul_small(M0, uint64_t(q));
    out.M256 = M0;
    if (!(tot == M0)) throw std::runtime_error("hier total != M0");
    return out;
}

static std::map<int, u128> brute(std::vector<int> const& primes) {
    int64_t M = 1;
    int T = 1;
    for (int q : primes) { M *= q; T = std::lcm(T, ord2(q)); }
    std::map<int, u128> h;
    for (int64_t a = 0; a < M; ++a) {
        int nu = 0;
        for (int r = 0; r < T; ++r) {
            bool ok = true;
            for (int q : primes) {
                int64_t p = 1 % q;
                for (int k = 0; k < r % ord2(q); ++k) p = p * 2 % q;
                if (a % q == p) { ok = false; break; }
            }
            nu += ok;
        }
        h[nu] += 1;
    }
    return h;
}

int main(int argc, char** argv) {
    std::vector<int> primes;
    std::vector<std::vector<int>> groups(1);
    size_t agroups_start = 0;
    bool self_test = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--self-test") self_test = true;
        else if (a == "--affine") {}  // always on (kept for command-line compatibility)
        else if (a == "--affine-max") AFFINE_MAX = std::stoi(argv[++i]);
        else if (a == "--fuse") FUSE = std::stoi(argv[++i]);
        else if (a == "--threads") THREADS = std::stoi(argv[++i]);
        else if (a == "--split") SPLIT = std::stoi(argv[++i]);
        else if (a == "--groups") GROUPS_G = std::stoi(argv[++i]);
        else if (a == "--stream") STREAM = std::stoi(argv[++i]);
        else if (a == "/") groups.emplace_back();
        else if (a == "--hier") HIER = true;
        else if (a == "--raw-limit") RAW_LIMIT = std::stoull(argv[++i]);
        else if (a == "--sigma-limit") SIGMA_LIMIT = std::stoull(argv[++i]);
        else if (a == "|") { agroups_start = groups.size(); groups.emplace_back(); }
        else { primes.push_back(std::stoi(a)); groups.back().push_back(primes.back()); }
    }
    try {
        if (self_test) {
            const int f0 = FUSE, t0 = THREADS, a0 = AFFINE_MAX;
            for (int am : {0, 420})
                for (int f : {0, 1, 2, 4})
                    for (int t : {1, 2}) {
                        AFFINE_MAX = am; FUSE = f; THREADS = t;
                        for (auto const& P : std::vector<std::vector<int>>{
                                 {3, 5, 7, 13}, {7, 3, 31, 5}, {17, 3, 5}, {13, 7, 5, 3}, {43, 127, 3, 7}}) {
                            if (run(P, false).h != brute(P)) throw std::runtime_error("self-test failed");
                        }
                    }
            for (int am : {0, 420})
                for (int f : {0, 2})
                    for (int k : {1, 2, 3}) {
                        AFFINE_MAX = am; FUSE = f; THREADS = 2;
                        for (auto const& P : std::vector<std::vector<int>>{
                                 {43, 3, 5, 7}, {127, 43, 3}, {29, 3, 5, 7}, {13, 17, 7, 3}, {7, 31, 3, 5}}) {
                            std::vector<int> A(P.begin(), P.begin() + k), B(P.begin() + k, P.end());
                            if (split_run(A, B, false).h != brute(P)) throw std::runtime_error("split self-test failed");
                        }
                    }
            for (int am : {0, 420}) {
                AFFINE_MAX = am; FUSE = 2; THREADS = 2;
                std::vector<std::vector<int>> G1{{3, 7, 5, 13}, {29}, {73, 19}};
                std::vector<std::vector<int>> G2{{29, 43}, {3, 5}, {73}, {17, 7}};
                for (auto const* GG : {&G1, &G2}) {
                    std::vector<int> all;
                    for (auto const& G : *GG) all.insert(all.end(), G.begin(), G.end());
                    FUSE = 0;
                    auto ref = run(all, false).h;
                    FUSE = 2;
                    for (int gg : {12, 4})
                        try {
                            for (int sp : {0, 3}) {
                                STREAM = sp;
                                if (group_run(*GG, gg, false).h != ref) throw std::runtime_error("group self-test failed");
                            }
                            STREAM = 0;
                        } catch (std::runtime_error const& e) {
                            if (std::string(e.what()).find("private parts") == std::string::npos) throw;
                        }
                }
            }
            {
                AFFINE_MAX = 420; FUSE = 2; THREADS = 2;
                std::vector<std::vector<int>> B{{3, 7, 5, 13}, {73, 19}, {17}}, A{{29, 43}, {23}};
                std::vector<std::vector<int>> all = B;
                all.insert(all.end(), A.begin(), A.end());
                if (hier_run(B, A, false).h != group_run(all, 12, false).h) throw std::runtime_error("hier self-test failed");
                std::vector<std::vector<int>> B2{{3, 7}, {31, 11}}, A2{{29}, {8191, 53}};
                std::vector<std::vector<int>> all2 = B2;
                all2.insert(all2.end(), A2.begin(), A2.end());
                if (hier_run(B2, A2, false).h != group_run(all2, 12, false).h) throw std::runtime_error("hier self-test 2 failed");
            }
            {
                U256 x = mul128(~u128(0), ~u128(0));  // (2^128-1)^2 = 2^256 - 2^129 + 1
                if (!(x.lo == 1 && x.hi == ~u128(0) - 1)) throw std::runtime_error("mul128 self-test failed");
                U256 y; y.lo = 1;
                for (int i = 0; i < 40; ++i) y = mul_small(y, 10);
                if (u256s(y) != "1" + std::string(40, '0')) throw std::runtime_error("u256s self-test failed");
                U256 z = mul128(u128(12345678901234567ull) * 1000000000000ull, u128(98765432109876543ull) * 1000000000ull);
                if (u256s(z) != "1219326311370217861743636654061881000000000000000000000")
                    throw std::runtime_error("u256 product self-test failed: " + u256s(z));
            }
            FUSE = f0; THREADS = t0; AFFINE_MAX = a0;
            std::cerr << "dp_seed2 self-test passed\n";
            if (primes.empty()) return 0;
        }
        if (HIER) {
            std::vector<std::vector<int>> B(groups.begin(), groups.begin() + agroups_start),
                A(groups.begin() + agroups_start, groups.end());
            Hist r = hier_run(B, A, true);
            for (auto const& [nu, c] : r.h256) std::cout << nu << " " << u256s(c) << "\n";
            std::cerr << "entries = " << r.h.size() << "\n";
            return 0;
        }
        if (GROUPS_G > 0) {
            Hist r = group_run(groups, GROUPS_G, true);
            for (auto const& [nu, c] : r.h) std::cout << nu << " " << u128s(c) << "\n";
            std::cerr << "entries = " << r.h.size() << "\n";
            return 0;
        }
        Hist r = SPLIT > 0 ? split_run(std::vector<int>(primes.begin(), primes.begin() + SPLIT),
                                       std::vector<int>(primes.begin() + SPLIT, primes.end()), true)
                           : run(primes, true);
        for (auto const& [nu, c] : r.h) std::cout << nu << " " << u128s(c) << "\n";
        std::cerr << "entries = " << r.h.size() << "\n";
    } catch (std::exception const& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}

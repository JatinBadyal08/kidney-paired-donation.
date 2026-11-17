"""Kidney Paired Donation: graph building, cycle detection, cycle packing."""
from collections import deque
import random

ABO = {"O": {"O", "A", "B", "AB"}, "A": {"A", "AB"}, "B": {"B", "AB"}, "AB": {"AB"}}
TYPES = ["O", "A", "B", "AB"]

def can_donate(donor, patient):
    return patient in ABO[donor]

def random_pairs(n):
    """Random (donor, patient, tissue) rows; each pair's own donor cannot give to its own patient."""
    out = []
    while len(out) < n:
        d, p = random.choice(TYPES), random.choice(TYPES)
        if not can_donate(d, p):
            out.append((d, p, random.randint(30, 98)))
    return out

def weight(a, b):
    """Tissue match score of a swap edge: average of the two pairs' scores (0-100)."""
    return (a[3] + b[3]) // 2

def build_matrix(pairs, min_score=0):
    """pairs: (id, donor, patient, tissue). Edge i->j needs ABO fit AND tissue score >= min_score."""
    n = len(pairs)
    return [[int(i != j and can_donate(pairs[i][1], pairs[j][2]) and weight(pairs[i], pairs[j]) >= min_score)
             for j in range(n)] for i in range(n)]

def cycle_score(c, pairs):
    return sum(weight(pairs[c[k]], pairs[c[(k + 1) % len(c)]]) for k in range(len(c)))

def two_cycles_matrix(adj):
    """Algorithm 1: matrix intersection, adj[a][b] == 1 and adj[b][a] == 1."""
    n = len(adj)
    return [(a, b) for a in range(n) for b in range(a + 1, n) if adj[a][b] and adj[b][a]]

def find_cycles(adj, max_len=3):
    """Algorithm 2: constrained DFS with backtracking (2- and 3-cycles).
    Only nodes with index > start are visited, so each cycle is reported once."""
    n, cycles = len(adj), []
    def dfs(start, node, path):
        for nxt in range(n):
            if not adj[node][nxt]:
                continue
            if nxt == start and len(path) >= 2:
                cycles.append(tuple(path))
            elif nxt > start and nxt not in path and len(path) < max_len:
                path.append(nxt)
                dfs(start, nxt, path)
                path.pop()
    for s in range(n):
        dfs(s, s, [s])
    return cycles

def greedy_packing(cycles, pairs):
    """Algorithm 3: shortest cycles first (best tissue score breaks ties), skip reused pairs."""
    chosen, used = [], set()
    for c in sorted(cycles, key=lambda c: (len(c), -cycle_score(c, pairs))):
        if used.isdisjoint(c):
            chosen.append(c)
            used.update(c)
    return chosen

def shortest_cycle(adj, s, max_len=3):
    """Algorithm 4: BFS with a FIFO queue; first loop back to s is the shortest."""
    q = deque([(s, [s])])
    while q:
        node, path = q.popleft()
        for nxt in range(len(adj)):
            if not adj[node][nxt]:
                continue
            if nxt == s and len(path) >= 2:
                return path
            if nxt not in path and len(path) < max_len:
                q.append((nxt, path + [nxt]))
    return None

def optimal_packing(cycles, pairs, limit=40):
    """Exhaustive search maximising (transplants, total tissue score); compared against greedy."""
    if len(cycles) > limit:
        return None
    best = [[], (0, 0)]
    def solve(i, chosen, used, tx, sc):
        if (tx, sc) > tuple(best[1]):
            best[0], best[1] = chosen[:], (tx, sc)
        if i == len(cycles):
            return
        c = cycles[i]
        if used.isdisjoint(c):
            chosen.append(c)
            solve(i + 1, chosen, used | set(c), tx + len(c), sc + cycle_score(c, pairs))
            chosen.pop()
        solve(i + 1, chosen, used, tx, sc)
    solve(0, [], set(), 0, 0)
    return best[0]

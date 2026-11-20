import sqlite3
from flask import Flask, jsonify, request, render_template
import algorithms as alg

app = Flask(__name__)
DB = "kpd.db"

def db():
    con = sqlite3.connect(DB)
    con.execute("CREATE TABLE IF NOT EXISTS pairs(id INTEGER PRIMARY KEY AUTOINCREMENT, donor TEXT, patient TEXT, tissue INTEGER DEFAULT 70)")
    try:  # upgrade databases created before tissue scores existed
        con.execute("ALTER TABLE pairs ADD COLUMN tissue INTEGER DEFAULT 70")
    except sqlite3.OperationalError:
        pass
    con.execute("CREATE TABLE IF NOT EXISTS compatibility(donor_pair INTEGER, patient_pair INTEGER)")
    return con

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/api/pairs", methods=["POST"])
def add_pair():
    j = request.json
    d, p = j.get("donor"), j.get("patient")
    try:
        t = int(j.get("tissue"))
    except (TypeError, ValueError):
        t = -1
    if d not in alg.TYPES or p not in alg.TYPES:
        return jsonify(error="Choose a blood type from O, A, B or AB."), 400
    if not 0 <= t <= 100:
        return jsonify(error="Tissue match score must be a number from 0 to 100."), 400
    if alg.can_donate(d, p):
        return jsonify(error=f"Donor {d} can already give to patient {p}, so this pair does not need a swap."), 400
    with db() as con:
        con.execute("INSERT INTO pairs(donor, patient, tissue) VALUES(?,?,?)", (d, p, t))
    return jsonify(ok=True)

@app.route("/api/pairs/<int:pid>", methods=["DELETE"])
def delete_pair(pid):
    with db() as con:
        con.execute("DELETE FROM pairs WHERE id=?", (pid,))
    return jsonify(ok=True)

@app.route("/api/generate", methods=["POST"])
def generate():
    n = max(2, min(int(request.json.get("n", 8)), 14))
    m = int(request.json.get("min", 0))
    for _ in range(300):  # retry until the random set contains at least one swap
        pairs = alg.random_pairs(n)
        rows = [(i, d, p, t) for i, (d, p, t) in enumerate(pairs)]
        if alg.find_cycles(alg.build_matrix(rows, m)):
            break
    with db() as con:
        con.execute("DELETE FROM pairs")
        con.executemany("INSERT INTO pairs(donor, patient, tissue) VALUES(?,?,?)", pairs)
    return jsonify(ok=True)

@app.route("/api/run")
def run():
    m = int(request.args.get("min", 0))
    with db() as con:
        pairs = con.execute("SELECT id, donor, patient, tissue FROM pairs ORDER BY id").fetchall()
        adj = alg.build_matrix(pairs, m)
        edges = [(pairs[i][0], pairs[j][0], alg.weight(pairs[i], pairs[j]))
                 for i in range(len(adj)) for j in range(len(adj)) if adj[i][j]]
        con.execute("DELETE FROM compatibility")
        con.executemany("INSERT INTO compatibility VALUES(?,?)", [e[:2] for e in edges])
    ids = lambda cs: [[pairs[i][0] for i in c] for c in cs]
    cycles = alg.find_cycles(adj)
    optimal = alg.optimal_packing(cycles, pairs)
    short = {pairs[i][0]: [pairs[k][0] for k in c] if (c := alg.shortest_cycle(adj, i)) else None
             for i in range(len(pairs))}
    return jsonify(
        pairs=[dict(id=i, donor=d, patient=p, tissue=t) for i, d, p, t in pairs],
        edges=edges,
        two_cycles=ids(alg.two_cycles_matrix(adj)),
        cycles=ids(cycles),
        greedy=ids(alg.greedy_packing(cycles, pairs)),
        optimal=None if optimal is None else ids(optimal),
        shortest=short)

if __name__ == "__main__":
    app.run()

import argparse
import csv
import json
import os
import sys

import iv_model

VALUE_KEYS = ["isc_stc_a", "i0_stc_a", "n_factor", "rs_ohm", "rsh_ohm",
              "ns_cells", "alpha_rel", "eg_ev", "noct_c", "t_amb_c", "g_wm2"]


def parse_row(row):
    if len(row) != 1 + len(VALUE_KEYS):
        raise ValueError("列数不对")
    b = {"batch_id": row[0]}
    b["isc_stc"] = float(row[1])
    b["i0"] = float(row[2])
    b["n"] = float(row[3])
    b["rs"] = float(row[4])
    b["rsh"] = float(row[5])
    b["ns"] = int(row[6])
    b["alpha"] = float(row[7])
    b["eg"] = float(row[8])
    b["noct"] = float(row[9])
    b["t_amb"] = float(row[10])
    b["g"] = float(row[11])
    return b


def model_params(batch):
    t_cell = iv_model.cell_temperature(batch["t_amb"], batch["noct"], batch["g"])
    nvt = iv_model.thermal_voltage(batch["n"], t_cell, batch["ns"])
    iph = iv_model.photocurrent(batch["isc_stc"], batch["alpha"], t_cell, batch["g"])
    i0 = iv_model.saturation_current(batch["i0"], batch["eg"], batch["n"], t_cell)
    return iph, i0, nvt


def evaluate(batch):
    iph, i0, nvt = model_params(batch)
    isc = iv_model.short_circuit_current(iph, i0, batch["rs"], batch["rsh"], nvt)
    voc = iv_model.open_circuit_voltage(iph, i0, batch["rs"], batch["rsh"], nvt)
    vmpp, pmpp = iv_model.max_power_point(iph, i0, batch["rs"], batch["rsh"], nvt)
    impp = pmpp / vmpp if vmpp > 0 else 0.0
    return isc, voc, vmpp, impp, pmpp


def write_curve(curve_dir, batch):
    iph, i0, nvt = model_params(batch)
    path = os.path.join(curve_dir, "%s.csv" % batch["batch_id"])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["v", "i"])
        for k in range(121):
            v = k * 0.5
            i = iv_model.current_at_voltage(v, iph, i0, batch["rs"], batch["rsh"], nvt)
            w.writerow(["%.1f" % v, "%.4f" % i])


def load_ledger(path):
    """读台账（按 batch_id 索引）；某行不是合法 JSON 时提示并把台账视为空，后面重建。"""
    records = {}
    if not os.path.exists(path):
        return records
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                records[rec["batch_id"]] = rec
            except (json.JSONDecodeError, KeyError, TypeError):
                print("台账 %s 有损坏行，本次已重建。" % path, file=sys.stderr)
                return {}
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--curve-dir", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--measurements")
    ap.add_argument("--fit-out")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
    os.makedirs(args.curve_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.ledger)), exist_ok=True)

    with open(args.input, newline="", encoding="utf-8") as f:
        rows = [row for row in csv.reader(f) if row][1:]

    results = {}
    order = []
    for row in rows:
        bid = row[0]
        order.append(bid)
        if bid in results:
            continue
        try:
            results[bid] = evaluate(parse_row(row))
        except (ValueError, TypeError):
            results[bid] = None

    with open(args.report, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["batch_id", "isc_a", "voc_v", "vmpp_v", "impp_a", "pmpp_w"])
        for bid in order:
            res = results[bid]
            if res is None:
                w.writerow([bid, "ERROR"])
            else:
                isc, voc, vmpp, impp, pmpp = res
                w.writerow([bid, "%.3f" % isc, "%.3f" % voc,
                            "%.3f" % vmpp, "%.3f" % impp, "%.3f" % pmpp])

    for row in rows:
        bid = row[0]
        if results.get(bid) is not None:
            try:
                write_curve(args.curve_dir, parse_row(row))
            except (ValueError, TypeError):
                pass

    if args.measurements and args.fit_out:
        import iv_fit
        cards = {}
        for row in rows:
            try:
                b = parse_row(row)
                cards[b["batch_id"]] = b
            except (ValueError, TypeError):
                continue
        with open(args.fit_out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["batch_id", "i0_a", "n_factor", "rs_ohm", "rsh_ohm", "rmse_a"])
            for bid, b in cards.items():
                mpath = os.path.join(args.measurements, "%s.csv" % bid)
                if not os.path.exists(mpath):
                    continue
                pts = []
                with open(mpath, newline="", encoding="utf-8") as mf:
                    for rrow in csv.DictReader(mf):
                        pts.append((float(rrow["v"]), float(rrow["i"])))
                t_cell = iv_model.cell_temperature(b["t_amb"], b["noct"], b["g"])
                iph = iv_model.photocurrent(b["isc_stc"], b["alpha"], t_cell, b["g"])
                i0, n, rs, rsh, rmse = iv_fit.identify(pts, iph, t_cell, b["ns"])
                w.writerow([bid, "%.3e" % i0, "%.4f" % n,
                            "%.4f" % rs, "%.4f" % rsh, "%.4f" % rmse])

    # 台账按 batch_id upsert：参数没变内容不变、不新增行；变了更新该行
    ledger = load_ledger(args.ledger)
    for bid in order:
        res = results[bid]
        if res is None:
            continue
        isc, voc, vmpp, impp, pmpp = res
        ledger[bid] = {"batch_id": bid, "isc_a": round(isc, 3),
                       "voc_v": round(voc, 3), "pmpp_w": round(pmpp, 3)}
    with open(args.ledger, "w", encoding="utf-8") as f:
        for bid in ledger:
            f.write(json.dumps(ledger[bid]) + "\n")


if __name__ == "__main__":
    main()

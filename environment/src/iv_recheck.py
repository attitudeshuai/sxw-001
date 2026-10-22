import argparse
import csv
import json
import os

import iv_model

VALUE_KEYS = ["isc_stc_a", "i0_stc_a", "n_factor", "rs_ohm", "rsh_ohm",
              "ns_cells", "alpha_rel", "eg_ev", "noct_c", "t_amb_c", "g_wm2"]


def read_batches(path):
    batches = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            batches.append({
                "batch_id": row["batch_id"],
                "isc_stc": float(row["isc_stc_a"]),
                "i0": float(row["i0_stc_a"]),
                "n": float(row["n_factor"]),
                "rs": float(row["rs_ohm"]),
                "rsh": float(row["rsh_ohm"]),
                "ns": int(row["ns_cells"]),
                "alpha": float(row["alpha_rel"]),
                "eg": float(row["eg_ev"]),
                "noct": float(row["noct_c"]),
                "t_amb": float(row["t_amb_c"]),
                "g": float(row["g_wm2"]),
            })
    return batches


def evaluate(batch):
    t_cell = iv_model.cell_temperature(batch["t_amb"], batch["noct"], batch["g"])
    nvt = iv_model.thermal_voltage(batch["n"], t_cell, batch["ns"])
    iph = iv_model.photocurrent(batch["isc_stc"], batch["alpha"], t_cell, batch["g"])
    i0 = iv_model.saturation_current(batch["i0"], batch["eg"], batch["n"], t_cell)
    isc = iv_model.short_circuit_current(iph, i0, batch["rs"], batch["rsh"], nvt)
    voc = iv_model.open_circuit_voltage(iph, i0, batch["rs"], batch["rsh"], nvt)
    vmpp, pmpp = iv_model.max_power_point(iph, i0, batch["rs"], batch["rsh"], nvt)
    impp = pmpp / vmpp
    return t_cell, iph, i0, isc, voc, vmpp, impp, pmpp


def write_curve(curve_dir, batch):
    t_cell = iv_model.cell_temperature(batch["t_amb"], batch["noct"], batch["g"])
    nvt = iv_model.thermal_voltage(batch["n"], t_cell, batch["ns"])
    iph = iv_model.photocurrent(batch["isc_stc"], batch["alpha"], t_cell, batch["g"])
    i0 = iv_model.saturation_current(batch["i0"], batch["eg"], batch["n"], t_cell)
    path = os.path.join(curve_dir, "%s.csv" % batch["batch_id"])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["v", "i"])
        for k in range(121):
            v = k * 0.5
            i = iv_model.current_at_voltage(v, iph, i0, batch["rs"], batch["rsh"], nvt)
            w.writerow(["%.1f" % v, "%.4f" % i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--curve-dir", required=True)
    ap.add_argument("--ledger", required=True)
    ap.add_argument("--measurements")
    ap.add_argument("--fit-out")
    args = ap.parse_args()

    batches = read_batches(args.input)

    os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
    os.makedirs(args.curve_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.ledger)), exist_ok=True)

    with open(args.report, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["batch_id", "isc_a", "voc_v", "vmpp_v", "impp_a", "pmpp_w"])
        for b in batches:
            _t, _iph, _i0, isc, voc, vmpp, impp, pmpp = evaluate(b)
            w.writerow([b["batch_id"], "%.3f" % isc, "%.3f" % voc,
                        "%.3f" % vmpp, "%.3f" % impp, "%.3f" % pmpp])

    for b in batches:
        write_curve(args.curve_dir, b)

    if args.measurements and args.fit_out:
        import iv_fit
        with open(args.fit_out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["batch_id", "i0_a", "n_factor", "rs_ohm", "rsh_ohm", "rmse_a"])
            for b in batches:
                mpath = os.path.join(args.measurements, "%s.csv" % b["batch_id"])
                if not os.path.exists(mpath):
                    continue
                pts = []
                with open(mpath, newline="", encoding="utf-8") as mf:
                    for row in csv.DictReader(mf):
                        pts.append((float(row["v"]), float(row["i"])))
                t_cell = iv_model.cell_temperature(b["t_amb"], b["noct"], b["g"])
                iph = iv_model.photocurrent(b["isc_stc"], b["alpha"], t_cell, b["g"])
                i0, n, rs, rsh, rmse = iv_fit.identify(pts, iph, t_cell, b["ns"])
                w.writerow([b["batch_id"], "%.3e" % i0, "%.4f" % n,
                            "%.4f" % rs, "%.4f" % rsh, "%.4f" % rmse])

    with open(args.ledger, "a", encoding="utf-8") as f:
        for b in batches:
            _t, _iph, _i0, isc, voc, vmpp, impp, pmpp = evaluate(b)
            f.write(json.dumps({"batch_id": b["batch_id"], "isc_a": round(isc, 3),
                                "voc_v": round(voc, 3), "pmpp_w": round(pmpp, 3)}) + "\n")


if __name__ == "__main__":
    main()

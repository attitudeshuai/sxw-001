import argparse
import csv
import json
import math
import os
import sys

import iv_model


def read_batches(path):
    batches = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            batches.append({
                "batch_id": row["batch_id"],
                "isc_stc": float(row["isc_stc_a"]),
                "i0": float(row["i0_a"]),
                "n": float(row["n_factor"]),
                "rs": float(row["rs_ohm"]),
                "rsh": float(row["rsh_ohm"]),
                "ns": int(row["ns_cells"]),
                "alpha": float(row["alpha_rel"]),
                "t_cell": float(row["t_cell_c"]),
                "g": float(row["g_wm2"]),
            })
    return batches


def evaluate(batch):
    nvt = iv_model.thermal_voltage(batch["n"], batch["t_cell"], batch["ns"])
    iph = iv_model.photocurrent(batch["isc_stc"], batch["alpha"], batch["t_cell"], batch["g"])
    isc = iv_model.short_circuit_current(iph, batch["i0"], batch["rs"], batch["rsh"], nvt)
    voc = iv_model.open_circuit_voltage(iph, batch["i0"], batch["rs"], batch["rsh"], nvt)
    vmpp, pmpp = iv_model.max_power_point(iph, batch["i0"], batch["rs"], batch["rsh"], nvt)
    impp = pmpp / vmpp
    return isc, voc, vmpp, impp, pmpp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--ledger", required=True)
    args = ap.parse_args()

    batches = read_batches(args.input)

    os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.ledger)), exist_ok=True)

    with open(args.report, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["batch_id", "isc_a", "voc_v", "vmpp_v", "impp_a", "pmpp_w"])
        for b in batches:
            isc, voc, vmpp, impp, pmpp = evaluate(b)
            w.writerow([b["batch_id"], "%.3f" % isc, "%.3f" % voc,
                        "%.3f" % vmpp, "%.3f" % impp, "%.3f" % pmpp])

    with open(args.ledger, "a", encoding="utf-8") as f:
        for b in batches:
            isc, voc, vmpp, impp, pmpp = evaluate(b)
            f.write(json.dumps({"batch_id": b["batch_id"], "isc_a": round(isc, 3),
                                "voc_v": round(voc, 3), "pmpp_w": round(pmpp, 3)}) + "\n")


if __name__ == "__main__":
    main()

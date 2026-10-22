import argparse
import csv
import json
import os
import sys

import iv_model

VALUE_KEYS = ["isc_stc_a", "i0_a", "n_factor", "rs_ohm", "rsh_ohm",
              "ns_cells", "alpha_rel", "t_cell_c", "g_wm2"]


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
    b["t_cell"] = float(row[8])
    b["g"] = float(row[9])
    return b


def evaluate(batch):
    nvt = iv_model.thermal_voltage(batch["n"], batch["t_cell"], batch["ns"])
    iph = iv_model.photocurrent(batch["isc_stc"], batch["alpha"], batch["t_cell"], batch["g"])
    isc = iv_model.short_circuit_current(iph, batch["i0"], batch["rs"], batch["rsh"], nvt)
    voc = iv_model.open_circuit_voltage(iph, batch["i0"], batch["rs"], batch["rsh"], nvt)
    vmpp, pmpp = iv_model.max_power_point(iph, batch["i0"], batch["rs"], batch["rsh"], nvt)
    impp = pmpp / vmpp if vmpp > 0 else 0.0
    return isc, voc, vmpp, impp, pmpp


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
    ap.add_argument("--ledger", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
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

    # 台账按 batch_id  upsert：参数没变时内容不变、不会新增行；变了就更新该批次的行
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

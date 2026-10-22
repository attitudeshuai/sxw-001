#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""判分：独立实现单二极管模型作为基准，驱动被测程序的 CLI 合同逐项核对。

按重要程度加权（不平均分）：核心数值 > 温度/辐照折算链 > 边界与错误处理 > 输出格式与曲线 > 台账。
基准侧的隐式方程用牛顿迭代（与中位二分不是同一条代码路径），曲线逐点核对。
"""
import csv
import json
import math
import os
import re
import subprocess
import sys
import tempfile

K_B = 1.380649e-23
Q_E = 1.602176634e-19
T_STC = 298.15

WEIGHTS = {
    "test_stc_case_matches_reference": 0.10,
    "test_analytic_voc_ideal_case": 0.10,
    "test_mpp_power_internally_consistent": 0.05,
    "test_voc_drops_when_cell_hot": 0.075,
    "test_isc_scales_with_irradiance": 0.075,
    "test_alpha_shifts_isc_off_stc": 0.05,
    "test_i0_temperature_shift_correct": 0.10,
    "test_noct_cell_temperature_used": 0.05,
    "test_cold_bright_batch_converges": 0.05,
    "test_huge_shunt_resistance_stable": 0.05,
    "test_illegal_row_marked_error": 0.05,
    "test_zero_irradiance_no_nan": 0.05,
    "test_report_header_and_columns": 0.05,
    "test_curve_full_grid_pointwise": 0.10,
    "test_ledger_rerun_no_duplicate": 0.02,
    "test_ledger_update_and_recover": 0.03,
}

HEADER = ["batch_id", "isc_a", "voc_v", "vmpp_v", "impp_a", "pmpp_w"]
BATCH_KEYS = ["isc_stc_a", "i0_stc_a", "n_factor", "rs_ohm", "rsh_ohm",
              "ns_cells", "alpha_rel", "eg_ev", "noct_c", "t_amb_c", "g_wm2"]

BASE = dict(isc_stc_a=9.20, i0_stc_a=1.0e-7, n_factor=1.20, rs_ohm=0.30, rsh_ohm=300.0,
            ns_cells=60, alpha_rel=0.00045, eg_ev=1.12, noct_c=45.0, t_amb_c=20.0, g_wm2=1000.0)

TOL_AV = 0.002
TOL_W = 0.02
TOL_CURVE = 0.005


# ---------------------------------------------------------------- 独立基准实现
def bench_params(b):
    t_cell = b["t_amb_c"] + (b["noct_c"] - 20.0) * b["g_wm2"] / 800.0
    t_k = t_cell + 273.15
    nvt = b["n_factor"] * K_B * t_k / Q_E * b["ns_cells"]
    iph = b["isc_stc_a"] * (1.0 + b["alpha_rel"] * (t_cell - 25.0)) * b["g_wm2"] / 1000.0
    i0 = b["i0_stc_a"] * (t_k / T_STC) ** 3 * math.exp(
        Q_E * b["eg_ev"] / (b["n_factor"] * K_B) * (1.0 / T_STC - 1.0 / t_k))
    return iph, i0, nvt


def i_at_v(v, iph, i0, rs, rsh, nvt):
    if iph <= 0:
        return 0.0
    # 牛顿迭代，近开路区切线步长过大时退回二分区间
    i = iph * 0.5
    lo, hi = -1e9, 1e9
    for _ in range(100):
        x = (v + i * rs) / nvt
        e = math.exp(x) if x < 700.0 else float("inf")
        f = iph - i0 * (e - 1.0) - (v + i * rs) / rsh - i
        df = -i0 * e * rs / nvt - rs / rsh - 1.0
        step = f / df if math.isfinite(df) and df != 0 else 0.0
        if f > 0:
            lo = i
        else:
            hi = i
        cand = i - step
        i = cand if lo < cand < hi else 0.5 * (lo + hi)
    return i


def voc_of(iph, i0, rs, rsh, nvt):
    if iph <= 0:
        return 0.0
    lo, hi = 0.0, 500.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if i_at_v(mid, iph, i0, rs, rsh, nvt) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def expected_points(b):
    iph, i0, nvt = bench_params(b)
    isc = i_at_v(0.0, iph, i0, b["rs_ohm"], b["rsh_ohm"], nvt)
    voc = voc_of(iph, i0, b["rs_ohm"], b["rsh_ohm"], nvt)
    best_v, best_p = 0.0, 0.0
    if voc > 0:
        step = voc / 400.0
        for k in range(401):
            v = voc * k / 400.0
            p = v * i_at_v(v, iph, i0, b["rs_ohm"], b["rsh_ohm"], nvt)
            if p > best_p:
                best_v, best_p = v, p
        v0, v2 = best_v - step, best_v + step
        if v0 >= 0.0 and v2 <= voc:
            p0 = v0 * i_at_v(v0, iph, i0, b["rs_ohm"], b["rsh_ohm"], nvt)
            p2 = v2 * i_at_v(v2, iph, i0, b["rs_ohm"], b["rsh_ohm"], nvt)
            denom = 2.0 * (p0 - 2.0 * best_p + p2)
            if abs(denom) > 1e-12:
                vref = best_v + (p0 - p2) * step / denom
                if v0 <= vref <= v2:
                    best_v = vref
                    best_p = best_v * i_at_v(best_v, iph, i0, b["rs_ohm"], b["rsh_ohm"], nvt)
    impp = best_p / best_v if best_v > 0 else 0.0
    return dict(isc=isc, voc=voc, vmpp=best_v, impp=impp, pmpp=best_p)


# ---------------------------------------------------------------- 驱动被测程序
SRC = os.environ.get("APP_SRC", "/app/src")


def run_tool(tmp, batches):
    inp = os.path.join(tmp, "in.csv")
    rep = os.path.join(tmp, "report.csv")
    led = os.path.join(tmp, "ledger.jsonl")
    cdir = os.path.join(tmp, "curve")
    with open(inp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["batch_id"] + BATCH_KEYS)
        for bid, b in batches:
            if isinstance(b, str):
                w.writerow([bid, b])   # 非法行
            else:
                w.writerow([bid] + [b[k] for k in BATCH_KEYS])
    r = subprocess.run([sys.executable, os.path.join(SRC, "iv_recheck.py"),
                        "--input", inp, "--report", rep, "--curve-dir", cdir, "--ledger", led],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=120)
    rows = []
    if os.path.exists(rep):
        with open(rep, newline="", encoding="utf-8") as f:
            rows = list(csv.reader(f))
    ledger = []
    if os.path.exists(led):
        with open(led, encoding="utf-8") as f:
            ledger = f.read().splitlines()
    return r, rows, ledger, cdir


def values(rows, bid):
    for row in rows[1:]:
        if row and row[0] == bid:
            if len(row) < 6 or row[1] == "ERROR":
                return None
            return dict(zip(("isc", "voc", "vmpp", "impp", "pmpp"), [float(x) for x in row[1:6]]))
    return None


def close(a, b, tol):
    return abs(a - b) <= tol


# ---------------------------------------------------------------- 测试
def test_stc_case_matches_reference():
    """题面对应原文：折算到现场工况，算出每块组件的关键点（短路电流、开路电压、最大功率点）。"""
    b = dict(BASE)
    exp = expected_points(b)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-STC", b)])
        got = values(rows, "T-STC")
    return (r.returncode == 0 and got
            and close(got["isc"], exp["isc"], TOL_AV) and close(got["voc"], exp["voc"], TOL_AV)
            and close(got["pmpp"], exp["pmpp"], TOL_W))


def test_analytic_voc_ideal_case():
    """题面对应原文：开路电压 Voc。"""
    b = dict(BASE, rs_ohm=0.0, rsh_ohm=1.0e12)
    iph, i0, nvt = bench_params(b)
    voc_ana = nvt * math.log(iph / i0 + 1.0)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-ANA", b)])
        got = values(rows, "T-ANA")
    return r.returncode == 0 and got and close(got["voc"], voc_ana, TOL_AV)


def test_mpp_power_internally_consistent():
    """题面对应原文：最大功率点（Vmpp、Impp、Pmpp）。"""
    b = dict(BASE, t_amb_c=30.0, g_wm2=900.0)
    exp = expected_points(b)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-MPP", b)])
        got = values(rows, "T-MPP")
    if not (r.returncode == 0 and got):
        return False
    return (close(got["vmpp"], exp["vmpp"], TOL_AV) and close(got["pmpp"], exp["pmpp"], TOL_W)
            and close(got["vmpp"] * got["impp"], got["pmpp"], 0.05))


def test_voc_drops_when_cell_hot():
    """题面对应原文：高温批次的开路电压总是偏高。"""
    hot = dict(BASE, t_amb_c=45.0)
    exp_base = expected_points(dict(BASE))
    exp_hot = expected_points(hot)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-20", dict(BASE)), ("T-45", hot)])
        g20, g45 = values(rows, "T-20"), values(rows, "T-45")
    return (r.returncode == 0 and g20 and g45
            and close(g45["voc"], exp_hot["voc"], TOL_AV)
            and close(g20["voc"] - g45["voc"], exp_base["voc"] - exp_hot["voc"], TOL_AV))


def test_isc_scales_with_irradiance():
    """题面对应原文：低辐照批次短路电流对不上手持表的读数。"""
    dim = dict(BASE, g_wm2=300.0)
    exp = expected_points(dim)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-G3", dim)])
        got = values(rows, "T-G3")
    return r.returncode == 0 and got and close(got["isc"], exp["isc"], TOL_AV)


def test_alpha_shifts_isc_off_stc():
    """题面对应原文：光电流 Iph = Isc_stc × (G/1000) × (1 + α×(T_cell − 25))。"""
    warm = dict(BASE, t_amb_c=35.0)
    exp = expected_points(warm)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-A35", warm)])
        got = values(rows, "T-A35")
    return r.returncode == 0 and got and close(got["isc"], exp["isc"], TOL_AV)


def test_i0_temperature_shift_correct():
    """题面对应原文：暗饱和电流随温度折算 I0(T) = I0_stc × (T/298.15)³ × exp[q·Eg/(n·k) × (1/298.15 − 1/T)]。"""
    cold = dict(BASE, t_amb_c=0.0, g_wm2=1000.0)
    exp = expected_points(cold)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-I0", cold)])
        got = values(rows, "T-I0")
    return r.returncode == 0 and got and close(got["voc"], exp["voc"], TOL_AV)


def test_noct_cell_temperature_used():
    """题面对应原文：电池温度按 NOCT 口径折算 T_cell = T_amb + (NOCT − 20) × G / 800。"""
    b = dict(BASE, t_amb_c=25.0, noct_c=48.0, g_wm2=800.0)
    exp = expected_points(b)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-NOCT", b)])
        got = values(rows, "T-NOCT")
    return r.returncode == 0 and got and close(got["isc"], exp["isc"], TOL_AV) and close(got["voc"], exp["voc"], TOL_AV)


def test_cold_bright_batch_converges():
    """题面对应原文：近开路区、零辐照、超大并联电阻等极端参数也要正常出数，不许死循环、溢出或 NaN。"""
    b = dict(BASE, isc_stc_a=9.60, i0_stc_a=8.0e-8, ns_cells=66, t_amb_c=-10.0, g_wm2=1100.0)
    exp = expected_points(b)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-CB", b)])
        got = values(rows, "T-CB")
    return (r.returncode == 0 and got
            and close(got["voc"], exp["voc"], TOL_AV) and close(got["pmpp"], exp["pmpp"], TOL_W))


def test_huge_shunt_resistance_stable():
    """题面对应原文：超大并联电阻等极端参数也要正常出数。"""
    b = dict(BASE, rsh_ohm=1.0e9)
    exp = expected_points(b)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-RSH", b)])
        got = values(rows, "T-RSH")
    return r.returncode == 0 and got and close(got["isc"], exp["isc"], TOL_AV)


def test_illegal_row_marked_error():
    """题面对应原文：非法行在报告里标 batch_id,ERROR，其余批次照常出数，退出码 0。"""
    good = dict(BASE)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-OK", good), ("T-BAD", "not-a-number")])
    ok_row = values(rows, "T-OK")
    bad = [row for row in rows[1:] if row and row[0] == "T-BAD"]
    exp = expected_points(good)
    return (r.returncode == 0 and ok_row and close(ok_row["isc"], exp["isc"], TOL_AV)
            and len(bad) == 1 and len(bad[0]) >= 2 and bad[0][1] == "ERROR")


def test_zero_irradiance_no_nan():
    """题面对应原文：不许死循环、溢出或 NaN。"""
    b = dict(BASE, g_wm2=0.0)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-G0", b)])
        got = values(rows, "T-G0")
    return (r.returncode == 0 and got
            and all(math.isfinite(v) for v in got.values())
            and close(got["isc"], 0.0, TOL_AV) and close(got["pmpp"], 0.0, TOL_W))


def test_report_header_and_columns():
    """题面对应原文：表头与列序固定；数值保留 3 位小数；批次顺序与输入一致。"""
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, _c = run_tool(tmp, [("T-H1", dict(BASE)), ("T-H2", dict(BASE, t_amb_c=35.0))])
    if r.returncode != 0 or not rows or rows[0] != HEADER:
        return False
    body = [row for row in rows[1:] if row and row[1] != "ERROR"]
    if [row[0] for row in body] != ["T-H1", "T-H2"]:
        return False
    return all(re.match(r"^-?\d+\.\d{3}$", c) for row in body for c in row[1:6])


def test_curve_full_grid_pointwise():
    """题面对应原文：电压固定取 0 到 60 V、每 0.5 V 一个点（共 121 点），每行 v,i，电流保留 4 位小数。"""
    b = dict(BASE, t_amb_c=30.0, g_wm2=850.0)
    iph, i0, nvt = bench_params(b)
    with tempfile.TemporaryDirectory() as tmp:
        r, rows, _l, cdir = run_tool(tmp, [("T-CV", b)])
        path = os.path.join(cdir, "T-CV.csv")
        if r.returncode != 0 or not os.path.exists(path):
            return False
        with open(path, newline="", encoding="utf-8") as f:
            curve = list(csv.reader(f))
    if not curve or curve[0] != ["v", "i"] or len(curve) != 122:
        return False
    for k, row in enumerate(curve[1:]):
        v = k * 0.5
        if abs(float(row[0]) - v) > 1e-9:
            return False
        exp_i = i_at_v(v, iph, i0, b["rs_ohm"], b["rsh_ohm"], nvt)
        if abs(float(row[1]) - exp_i) > TOL_CURVE:
            return False
    return True


def test_ledger_rerun_no_duplicate():
    """题面对应原文：同一批次参数没变时重跑不得新增行。"""
    with tempfile.TemporaryDirectory() as tmp:
        run_tool(tmp, [("T-L1", dict(BASE)), ("T-L2", dict(BASE, t_amb_c=30.0))])
        _r, _rows, led, _c = run_tool(tmp, [("T-L1", dict(BASE)), ("T-L2", dict(BASE, t_amb_c=30.0))])
    ids = [json.loads(x)["batch_id"] for x in led if x.strip()]
    return len(ids) == 2 and sorted(ids) == ["T-L1", "T-L2"]


def test_ledger_update_and_recover():
    """题面对应原文：参数变了更新该批次的行；台账有损坏行时给出可读提示并重建，退出码 0，不抛 traceback。"""
    with tempfile.TemporaryDirectory() as tmp:
        _r, _rows, _l, _c = run_tool(tmp, [("T-U1", dict(BASE))])
        led_path = os.path.join(tmp, "ledger.jsonl")
        changed = dict(BASE, t_amb_c=40.0)
        r2, _rows2, led2, _c2 = run_tool(tmp, [("T-U1", changed)])
        exp = expected_points(changed)
        upd_ok = False
        if r2.returncode == 0 and len([x for x in led2 if x.strip()]) == 1:
            rec = json.loads(led2[0])
            upd_ok = close(rec["voc_v"], round(exp["voc"], 3), TOL_AV)
        with open(led_path, "a", encoding="utf-8") as f:
            f.write("{broken json line\n")
        r3 = subprocess.run([sys.executable, os.path.join(SRC, "iv_recheck.py"),
                             "--input", os.path.join(tmp, "in.csv"),
                             "--report", os.path.join(tmp, "r3.csv"),
                             "--curve-dir", os.path.join(tmp, "c3"),
                             "--ledger", led_path],
                            capture_output=True, text=True, encoding="utf-8", errors="replace",
                            timeout=120)
        recovered = False
        if r3.returncode == 0 and "Traceback" not in r3.stderr:
            with open(led_path, encoding="utf-8") as f:
                lines = [x for x in f.read().splitlines() if x.strip()]
            try:
                ids = [json.loads(x)["batch_id"] for x in lines]
                recovered = ids == ["T-U1"] and bool((r3.stdout + r3.stderr).strip())
            except json.JSONDecodeError:
                recovered = False
    return upd_ok and recovered


def main():
    total = 0.0
    for name, w in WEIGHTS.items():
        fn = globals().get(name)
        try:
            passed = bool(fn and fn())
        except Exception as e:  # noqa: BLE001
            print("%s: 异常 %r" % (name, e), file=sys.stderr)
            passed = False
        total += w * passed
        print("%s: %s (权重 %.2f)" % (name, "PASS" if passed else "FAIL", w), file=sys.stderr)
    print("%.6f" % total)


if __name__ == "__main__":
    main()

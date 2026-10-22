"""单二极管模型的 IV 特性计算（修正版）。"""

import math

K_B = 1.380649e-23
Q_E = 1.602176634e-19


def thermal_voltage(n_factor, t_cell_c, ns_cells):
    # 热电压随电池温度走：Vt = n·k·T/q（T 用开尔文），组件再乘串联片数
    return n_factor * K_B * (t_cell_c + 273.15) / Q_E * ns_cells


def photocurrent(isc_stc, alpha_rel, t_cell_c, g_wm2):
    # 光电流随辐照线性缩放，温度项是相对 25°C 的增量
    return isc_stc * (1.0 + alpha_rel * (t_cell_c - 25.0)) * g_wm2 / 1000.0


def current_at_voltage(v, iph, i0, rs, rsh, nvt):
    if iph <= 0:
        return 0.0
    # 隐式方程对 I 单调，用二分法定摄；不动点迭代在近开路区会发散

    def gap(i):
        x = (v + i * rs) / nvt
        e = math.exp(x) if x < 700.0 else float("inf")
        return iph - i0 * (e - 1.0) - (v + i * rs) / rsh - i

    hi = iph  # gap(hi) <= 0
    lo, step = iph - 1.0, 1.0
    while gap(lo) <= 0:  # 向下扩括号；指数项随 i 下降快速衰减，几步就到位
        lo -= step
        step *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if gap(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def open_circuit_voltage(iph, i0, rs, rsh, nvt):
    if iph <= 0:
        return 0.0
    lo, hi = 0.0, 500.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if current_at_voltage(mid, iph, i0, rs, rsh, nvt) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def short_circuit_current(iph, i0, rs, rsh, nvt):
    return current_at_voltage(0.0, iph, i0, rs, rsh, nvt)


def max_power_point(iph, i0, rs, rsh, nvt):
    voc = open_circuit_voltage(iph, i0, rs, rsh, nvt)
    if voc <= 0:
        return 0.0, 0.0
    # P(V) 单峰，用黄金分割搜索替代粗扫
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = 0.0, voc

    def power(v):
        return v * current_at_voltage(v, iph, i0, rs, rsh, nvt)

    c = b - gr * (b - a)
    d = a + gr * (b - a)
    while abs(b - a) > 1e-7:
        if power(c) < power(d):
            a = c
        else:
            b = d
        c = b - gr * (b - a)
        d = a + gr * (b - a)
    v = 0.5 * (a + b)
    return v, power(v)

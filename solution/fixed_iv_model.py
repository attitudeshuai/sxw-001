"""单二极管模型的 IV 特性计算（修正版）。"""

import math

K_B = 1.380649e-23
Q_E = 1.602176634e-19
T_STC = 298.15


def cell_temperature(t_amb_c, noct_c, g_wm2):
    # NOCT 口径：基准工况是 800 W/m²、20°C
    return t_amb_c + (noct_c - 20.0) * g_wm2 / 800.0


def thermal_voltage(n_factor, t_cell_c, ns_cells):
    return n_factor * K_B * (t_cell_c + 273.15) / Q_E * ns_cells


def photocurrent(isc_stc, alpha_rel, t_cell_c, g_wm2):
    return isc_stc * (1.0 + alpha_rel * (t_cell_c - 25.0)) * g_wm2 / 1000.0


def saturation_current(i0_stc, eg_ev, n_factor, t_cell_c):
    # I0(T) = I0_stc·(T/Tstc)³·exp[q·Eg/(n·k)·(1/Tstc − 1/T)]，Eg 从 eV 换成焦耳
    t = t_cell_c + 273.15
    return i0_stc * (t / T_STC) ** 3 * math.exp(
        Q_E * eg_ev / (n_factor * K_B) * (1.0 / T_STC - 1.0 / t))


def current_at_voltage(v, iph, i0, rs, rsh, nvt):
    if iph <= 0:
        return 0.0

    def gap(i):
        x = (v + i * rs) / nvt
        e = math.exp(x) if x < 700.0 else float("inf")
        return iph - i0 * (e - 1.0) - (v + i * rs) / rsh - i

    hi = iph  # gap(hi) <= 0
    lo, step = iph - 1.0, 1.0
    while gap(lo) <= 0:  # 向下扩括号，指数项随 i 下降快速衰减
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

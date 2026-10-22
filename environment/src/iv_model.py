"""单二极管模型的 IV 特性计算。"""

import math

K_B = 1.380649e-23
Q_E = 1.602176634e-19
T_REF = 298.15


def thermal_voltage(n_factor, t_cell_c, ns_cells):
    v_t = n_factor * K_B * T_REF / Q_E
    return v_t * ns_cells


def photocurrent(isc_stc, alpha_rel, t_cell_c, g_wm2):
    return isc_stc * (1.0 + alpha_rel * t_cell_c) * g_wm2 / 1000.0


def current_at_voltage(v, iph, i0, rs, rsh, nvt):
    i = iph
    for _ in range(50):
        i = iph - i0 * (math.exp((v + i * rs) / nvt) - 1.0) - (v + i * rs) / rsh
    return i


def open_circuit_voltage(iph, i0, rs, rsh, nvt):
    lo, hi = 0.0, 200.0
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
    best_v, best_p = 0.0, 0.0
    v = 0.0
    while v <= voc:
        i = current_at_voltage(v, iph, i0, rs, rsh, nvt)
        p = v * i
        if p > best_p:
            best_v, best_p = v, p
        v += voc / 200.0
    return best_v, best_p

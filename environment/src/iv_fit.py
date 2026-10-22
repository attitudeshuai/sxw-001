"""参数辨识：从标定曲线估 (I0, n, Rs, Rsh)。"""

import math

import iv_model

K_B = 1.380649e-23
Q_E = 1.602176634e-19


def identify(points, iph, t_cell_c, ns_cells):
    n, rs, rsh, i0 = 1.3, 0.3, 300.0, 1e-6
    steps = {"n": 0.05, "rs": 0.05, "rsh": 20.0, "i0": 1e-7}

    def rmse(n_, rs_, rsh_, i0_):
        nvt = n_ * K_B * (t_cell_c + 273.15) / Q_E * ns_cells
        err = 0.0
        for v, i_meas in points:
            i_hat = iv_model.current_at_voltage(v, iph, i0_, rs_, rsh_, nvt)
            err += (i_hat - i_meas) ** 2
        return math.sqrt(err / len(points))

    cur = rmse(n, rs, rsh, i0)
    for _ in range(60):
        for name in ("n", "rs", "rsh", "i0"):
            val = {"n": n, "rs": rs, "rsh": rsh, "i0": i0}[name]
            for sign in (1, -1):
                cand = val + sign * steps[name]
                c = {"n": n, "rs": rs, "rsh": rsh, "i0": i0}
                c[name] = cand
                score = rmse(c["n"], c["rs"], c["rsh"], c["i0"])
                if score < cur:
                    cur = score
                    n, rs, rsh, i0 = c["n"], c["rs"], c["rsh"], c["i0"]
                    val = cand
    return i0, n, rs, rsh, cur

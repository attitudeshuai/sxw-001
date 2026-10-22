"""参数辨识（反问题）：从标定曲线反推 (I0, n, Rs, Rsh)。

做法：I0 与 Rsh 取 log10 拉到与其它参数同量级，Nelder-Mead 单纯形下山；
每个评价点上电流由隐式方程定摄（iv_model.current_at_voltage）。
"""

import math

import iv_model

K_B = 1.380649e-23
Q_E = 1.602176634e-19


def _unpack(x):
    log_i0, n, rs, log_rsh = x
    return 10.0 ** log_i0, n, rs, 10.0 ** log_rsh


def identify(points, iph, t_cell_c, ns_cells):
    """points: [(v, i)]；iph 为该工况的光电流（已知）；返回 (i0, n, rs, rsh, rmse)。"""
    t_k = t_cell_c + 273.15

    def objective(x):
        i0, n, rs, rsh = _unpack(x)
        if not (0.5 <= n <= 3.0) or rs < 0 or rs > 5 or rsh < 1 or rsh > 1e6:
            return 1e9
        nvt = n * K_B * t_k / Q_E * ns_cells
        err = 0.0
        for v, i_meas in points:
            i_hat = iv_model.current_at_voltage(v, iph, i0, rs, rsh, nvt)
            err += (i_hat - i_meas) ** 2
        return math.sqrt(err / len(points))

    x0 = [-6.0, 1.3, 0.3, 2.5]
    steps = [0.3, 0.1, 0.1, 0.3]
    simplex = [list(x0)]
    for j in range(4):
        v = list(x0)
        v[j] += steps[j]
        simplex.append(v)

    vals = [objective(x) for x in simplex]
    for _ in range(2000):
        order = sorted(range(5), key=lambda i: vals[i])
        simplex = [simplex[i] for i in order]
        vals = [vals[i] for i in order]
        spread = max(abs(v - vals[0]) for v in vals)
        size = max(max(abs(simplex[i][j] - simplex[0][j]) for j in range(4))
                   for i in range(1, 5))
        if spread < 1e-12 and size < 1e-8:
            break
        centroid = [sum(simplex[i][j] for i in range(4)) / 4.0 for j in range(4)]
        worst = simplex[4]
        xr = [centroid[j] + (centroid[j] - worst[j]) for j in range(4)]
        fr = objective(xr)
        if vals[0] <= fr < vals[3]:
            simplex[4], vals[4] = xr, fr
        elif fr < vals[0]:
            xe = [centroid[j] + 2.0 * (xr[j] - centroid[j]) for j in range(4)]
            fe = objective(xe)
            if fe < fr:
                simplex[4], vals[4] = xe, fe
            else:
                simplex[4], vals[4] = xr, fr
        else:
            xc = [centroid[j] + 0.5 * (worst[j] - centroid[j]) for j in range(4)]
            fc = objective(xc)
            if fc < vals[4]:
                simplex[4], vals[4] = xc, fc
            else:
                for i in range(1, 5):
                    simplex[i] = [simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j])
                                  for j in range(4)]
                    vals[i] = objective(simplex[i])

    best = simplex[vals.index(min(vals))]
    i0, n, rs, rsh = _unpack(best)
    return i0, n, rs, rsh, min(vals)

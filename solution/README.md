# 难点与测试对应

## 难点 1：温度折算链（NOCT、热电压、光电流、暗饱和电流）
推理轴：领域公式
落在推理路径上的哪里：现场工况折算是一整条链：T_cell 要按 NOCT 口径从环境温度折算（分母是 800 不是 1000）；热电压用电池温度的开尔文值；光电流的温度项基准是 25°C；暗饱和电流还要按 (T/298.15)³·exp[q·Eg/(n·k)·(1/298.15 − 1/T)] 走，Eg 要从 eV 换成焦耳、n 在指数分母里。原有实现四处全错：T_cell 折算分母错、热电压写死 25°C、温度项基准写成 0°C、I0 完全没折算。每一处独立影响不同工况的答案。
对应测试：test_noct_cell_temperature_used、test_voc_drops_when_cell_hot、test_alpha_shifts_isc_off_stc、test_i0_temperature_shift_correct

## 难点 2：隐式方程的全区间定摄
推理轴：数值方法
落在推理路径上的哪里：I = f(V, I) 隐式，近开路区不动点迭代振荡到溢出；而且合同要输出 0–60 V 固定网格的整条曲线（121 点、可为负值），关键点修好了、远区发散照样过不了逐点验收。要换成有 bracket 的方法并对指数溢出做保护。
对应测试：test_cold_bright_batch_converges、test_zero_irradiance_no_nan、test_huge_shunt_resistance_stable、test_curve_full_grid_pointwise

## 难点 3：台账与容错的工程合同
推理轴：状态迁移
落在推理路径上的哪里：台账按 batch_id upsert：重跑不新增、参数变了更新原行、损坏行提示后重建且退出码 0；非法行在报告里标 ERROR 但不拖垮其它批次。
对应测试：test_ledger_rerun_no_duplicate、test_ledger_update_and_recover、test_illegal_row_marked_error

## 跨运行状态
见难点 3：重复跑、参数变更、损坏恢复三种迁移都有独立测试（test_ledger_rerun_no_duplicate / test_ledger_update_and_recover 的更新与恢复两段）。

## 输出契约（另有独立测试）
报告表头列序与三位小数（test_report_header_and_columns）、曲线网格与四位小数（test_curve_full_grid_pointwise）、STC 参照与解析极限自检（test_stc_case_matches_reference、test_analytic_voc_ideal_case、test_mpp_power_internally_consistent）。

## 权重口径
核心数值 0.25、温度折算链 0.35、边界容错 0.20、格式与曲线 0.15、台账 0.05；权重在 tests/test_outputs.py 的 WEIGHTS 里，和为 1。

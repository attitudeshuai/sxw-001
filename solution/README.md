# 难点与测试对应

## 难点 1：热电压必须随电池温度走
推理轴：领域公式
落在推理路径上的哪里：单二极管模型里指数项的电压尺度是 n·k·T/q，T 要用电池温度的开尔文值；题面症状是高温批次开路电压偏高。原有实现把这个尺度写死在 25°C，温度一进参数就被忽略。
对应测试：test_voc_drops_when_cell_hot、test_alpha_shifts_isc_off_stc

## 难点 2：隐式方程的定摄方法
推理轴：数值方法
落在推理路径上的哪里：电流满足 I = f(V, I) 的隐式关系，近开路区不动点迭代来回振荡直至溢出；要换成对单调函数的二分（或等价的有 bracket 方法），并对指数溢出做保护。
对应测试：test_cold_bright_batch_converges、test_zero_irradiance_no_nan、test_huge_shunt_resistance_stable

## 难点 3：光电流的温度系数基准与辐照归一
推理轴：数据口径
落在推理路径上的哪里：光电流 = STC 短路电流 × 辐照/1000 × (1 + α·(T−25))；原有实现把温度项的基准从 25°C 弄成了 0°C，低辐照和高温批次全都对不上手持表。
对应测试：test_isc_scales_with_irradiance、test_alpha_shifts_isc_off_stc、test_stc_case_matches_reference

## 跨运行状态（另有两条独立守格式与容错）
台账按批次 upsert：重跑不重复登记、参数变了更新原行、文件损坏给可读提示并重建且退出码 0。对应测试：test_ledger_rerun_no_duplicate、test_ledger_update_and_recover。边界容错与输出契约另有 test_illegal_row_marked_error、test_report_header_and_columns、test_analytic_voc_ideal_case、test_mpp_power_internally_consistent 单独覆盖。

## 权重口径
核心数值 0.40（stc_case / analytic_voc / mpp_consistency），温度与辐照规则 0.25，边界与错误处理 0.20，输出格式 0.10，台账 0.05；权重在 tests/test_outputs.py 的 WEIGHTS 里，和为 1。

# 光伏组件 IV 曲线复核单

电站运维夜班要复核一批光伏组件的 IV 特性。厂家参数卡给的是单二极管模型参数（STC：25°C、1000 W/m² 标定）。班组的复核工具读批次清单，把参数折算到现场工况，算出每块组件的关键点（短路电流、开路电压、最大功率点），输出整条 IV 曲线，出复核报告，并登记台账。

## 折算口径（参数卡附带的工程规范，逐条照做）

- 电池温度按 NOCT 口径从环境温度折算：T_cell = T_amb + (NOCT − 20) × G / 800，NOCT 与 T_amb 单位 °C，G 单位 W/m²。
- 热电压用电池温度的开尔文值：Vt（每片）= k·T_cell/q，组件再乘串联片数与理想因子。
- 光电流：Iph = Isc_stc × (G/1000) × (1 + α×(T_cell − 25))，α 是相对温度系数。
- 暗饱和电流随温度折算：I0(T) = I0_stc × (T/298.15)³ × exp[ q·Eg/(n·k) × (1/298.15 − 1/T) ]，T 为 T_cell 开尔文值，Eg 单位 eV（公式里要换成焦耳），n 是理想因子。
- 单二极管方程：I = Iph − I0×(exp((V + I×Rs)/(n×Ns×Vt)) − 1) − (V + I×Rs)/Rsh。

## 现状

工具在 `src/` 下，跑法：

```
python3 src/iv_recheck.py --input data/batches.csv --report out/report.csv --curve-dir out/curve --ledger out/ledger.jsonl --measurements data/measurements --fit-out out/fit.csv
```

最近几班反映：高温批次的开路电压总是偏高；低辐照批次短路电流对不上手持表；近开路区程序会振荡甚至溢出中断；同一批次重复跑会重复登记台账。

## 任务

把复核工具改对，并让输出符合下面的合同。另外，有一批标定曲线要做参数辨识（见下）。

## 参数辨识（反问题）

`data/measurements/` 下有若干 `{batch_id}.csv` 标定曲线（每行 `v,i`，是运维平台按单二极管模型加固定生成规则产出的标定数据，覆盖 0 到 60 V）。这些批次的工况在 `data/batches.csv` 里查得到。对每个有标定曲线的批次，在该批次自己的工况下辨识四个参数：暗饱和电流 I0（A）、理想因子 n、串联电阻 Rs（Ω）、并联电阻 Rsh（Ω）。该工况下的光电流按折算口径先算出来，不是待辨识量。

结果写到 `--fit-out` 指定的 CSV，表头与列序固定：`batch_id,i0_a,n_factor,rs_ohm,rsh_ohm,rmse_a`。rmse_a 是用辨识出的参数回放曲线、与标定点电流的均方根误差（A）。辨识值保留到足够精度（I0 用科学计数法，其余 4 位小数）。没有标定曲线的批次不出现在 fit 结果里。

提醒：这四个参数的量级差得很远（I0 在 1e-7 量级，Rsh 到几百），直接在同一尺度上搜索是搜不动的。

## 输入

`data/batches.csv`，每行一个批次，表头顺序固定：

- batch_id：批次号
- isc_stc_a：STC 短路电流，A
- i0_stc_a：STC 暗饱和电流，A
- n_factor：二极管理想因子
- rs_ohm：串联电阻，Ω
- rsh_ohm：并联（旁路）电阻，Ω
- ns_cells：串联电池片数
- alpha_rel：短路电流相对温度系数，1/°C
- eg_ev：材料带隙，eV
- noct_c：额定工作温度 NOCT，°C
- t_amb_c：现场环境温度，°C
- g_wm2：现场辐照，W/m²

## 输出

1. 报告写到 `--report` 指定的 CSV，表头与列序固定：`batch_id,isc_a,voc_v,vmpp_v,impp_a,pmpp_w`，数值保留 3 位小数，批次顺序与输入一致。
2. 曲线写到 `--curve-dir` 目录下 `{batch_id}.csv`：电压固定取 0 到 60 V、每 0.5 V 一个点（共 121 点），每行 `v,i`，电流保留 4 位小数，允许负值。非法批次不产生曲线文件。
3. 台账写到 `--ledger` 指定的 JSONL，每行一个 JSON 对象，字段为 batch_id、isc_a、voc_v、pmpp_w：同一批次参数没变时重跑不得新增行；参数变了更新该批次的行；台账有损坏行时给出可读提示并重建，退出码 0，不抛 traceback。

## 约束

- 只用 Python 标准库。
- 非法行（列数不对、数值无法解析）在报告里标 `batch_id,ERROR`，其余批次照常出数，退出码 0。
- 近开路区、零辐照、超大并联电阻等极端参数也要正常出数，不许死循环、溢出或 NaN。

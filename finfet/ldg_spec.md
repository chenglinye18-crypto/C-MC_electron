# `ldg.txt` 配置说明（基于当前 FinFET 示例）

`ldg.txt` 是器件几何、材料、边界与接触的脚本式描述文件。程序按行解析关键字，顺序通常无严格要求，但某些命令相互依赖（如 `attachcontact` 与 `contact`）。以下按常见关键字说明用法与物理含义，并结合当前文件片段示例。

## 坐标与范围约定

- 每条命令都以 6 个坐标（`xbegin xend ybegin yend zbegin zend`）定义一个长方体或平面范围。坐标值需与 `lgrid.txt` 的网格坐标一致。
- 范围是闭区间，单位与 `lgrid.txt` 相同（原始多为微米，代码内乘 `1e-6` 变为米）。
- DIR 方向枚举：`UP(-X) DOWN(+X) LEFT(-Y) RIGHT(+Y) FRONT(-Z) BACK(+Z)`。

## 关键字列表

### default_par_number
```
default_par_number ele_min hole_min
```
设置每个 cell 的电子/空穴数下限，用于粒子初始化（未被其他命令覆盖时采用）。

示例：
```
default_par_number 15 4
```

### region
```
region x1 x2 y1 y2 z1 z2 MATERIAL
```
定义材料区域，同时初始化该区域的目标载流子数。常用 MATERIAL：`SILICON`、`OXIDE`、`VACUUM`。可多次覆盖。

示例：
```
region  0 60 -32 32 0 10 SILICON
region 60 160 -32 32 -30 40 SILICON
region -1 0 -32 32 -1 11 OXIDE
```

### donor / acceptor
```
donor    x1 x2 y1 y2 z1 z2 ND
acceptor x1 x2 y1 y2 z1 z2 NA
```
掺杂浓度分布（单位遵循程序内部常量，通常 m^-3）。`0` 表示无额外掺杂。

示例：
```
donor    0 40 -32 -10 0 10 1e20
acceptor 40 60 -32 32 0 10 2e16
```

### motioncube
```
motioncube x1 x2 y1 y2 z1 z2 MX+ MX- MY+ MY- MZ+ MZ-
```
为长方体六个面设置粒子运动规则（边界条件）。规则可为：`PASS`、`REFLECT`、`SCATTOX`、`CATCH`、`GENERATE`、`GENREF`、`PERIOD` 等。
- `PASS`：直接穿透到相邻单元。
- `REFLECT`：镜面反射。
- `SCATTOX`：Si/SiO2 界面散射（概率漫反/镜反）。
- `CATCH`：接触吸收。
- `GENERATE`/`GENREF`：在边界内侧生成一个拷贝粒子，同时当前粒子被反射（用于接触注入/周期补充）。
- `PERIOD`：周期边界。

示例（整体设为通过）：
```
motioncube 0 160 -32 32 -30 40 PASS PASS PASS PASS PASS PASS
```

### motionplane
```
motionplane x1 x2 y1 y2 z1 z2 DIR RULE
```
为单个平面设置方向敏感的运动规则。常用于源/漏接触与外层反射。

示例（左/右接触与外层反射）：
```
motionplane  0 40 -31 -31 0 10 LEFT  CATCH     # 粒子从 -Y 侧撞上被吸收
motionplane  0 40 -31 -31 0 10 RIGHT GENERATE  # 粒子从 +Y 侧撞上生成并反射
motionplane  0 60 -32 -32 0 10 LEFT  REFLECT   # 外层反射
```

### attachcontact
```
attachcontact x1 x2 y1 y2 z1 z2 INDEX
```
给硅区标注接触编号（INDEX>0）。后续 `contact` 按出现顺序对应 INDEX-1。用来让粒子知道落在哪个接触。

示例：
```
attachcontact 0 40 -32 32 0 10 0    # GATE=0
attachcontact 0 40 -32 -30 0 10 1   # SOURCE=1
attachcontact 0 40 30 32 0 10 2     # DRAIN=2
attachcontact 150 160 -32 32 -30 40 3 # BACK=3
```

### contact
块状定义接触的平面位置与施加电压/功函数差：
```
contact
  N phi_ms
  x1 x2 y1 y2 z1 z2   # 共 N 行
  ...
  Vapp                # 末行：接触电压
```
块的顺序决定索引：第 1 个 contact 对应 attachcontact 中的 1（内部索引 0），第 2 个对应 attachcontact 中的 2，依此类推。`phi_ms` 会除以 pot0 归一化。

示例（源/漏各 1 个平面，栅有 3 个平面）：
```
contact
  1 0
  0 40 -32 -32 0 10
  0           # Vapp

contact
  1 0
  0 40 32 32 0 10
  1           # Vapp

contact
  3 -0.46
  -1 -1 -10 10 -1 11
  -1 40 -10 10 -1 -1
  -1 40 -10 10 11 11
  1           # Vapp
```

### ScatterArea
```
ScatterArea x1 x2 y1 y2 z1 zend TYPE
```
标记散射区域及类型编号（具体含义在代码中自定义）。

示例：
```
ScatterArea 0 40 -32 32 0 10 3
```

### VsVdVg / vgrange
标记为 “not used anymore”，通常忽略。

### ep_parm
```
ep_parm qc_xratio qc_xtheta qc_yratio qc_ytheta qc_zratio qc_ztheta Eb
```
量子/界面参数：`qc_*theta` 为特征长度（nm 级，内部归一化），`qc_*ratio` 为乘子，`Eb` 为势垒高度（eV，内部除 pot0）。用于设置量子修正范围及氧化层势垒。

示例（默认）：
```
ep_parm 4 0.5 4 0.5 4 0.5 3.1
```

### parnumber
```
parnumber x1 x2 y1 y2 z1 z2 e_num h_num
```
为指定区域覆盖/设置目标电子、空穴数（替代 default_par_number 的默认值）。

示例：
```
parnumber 40 60 -32 32 0 10 5 5
```

### surfaces / surface_scatter_range
```
surfaces N
  surf_type surf_pos surf_dir  # 共 N 行
surface_scatter_range x1 x2 y1 y2 z1 z2
```
定义 Si/SiO2 界面类型与位置（`surf_type`: 0 yz, 1 xz, 2 xy），`surf_dir` 为界面法向指向氧化物。`surface_scatter_range` 指定表面散射作用范围。

示例：
```
surfaces 4
  0 0 0
  0 10 1
  2 0 4
  2 10 5
surface_scatter_range 0 40 -32 32 0 10
```

### quantumRegion
```
quantumRegion x1 x2 y1 y2 z1 z2
```
标记量子修正作用区域。

## 解析与注意事项
- 文件不支持注释行，避免在命令行后添加说明文字。
- `attachcontact` 的 INDEX 与 `contact` 块的顺序必须一致（内部使用 INDEX-1 存取）。
- 坐标要与 `lgrid.txt` 网格对应；平邻/重叠范围允许，但要与材料/接触物理一致。
- 规则选择建议：接触面常用 `CATCH`/`GENERATE` 成对；外层用 `REFLECT` 保护；Si/SiO2 界面用 `SCATTOX` 或 `REFLECT`；一般边界用 `PASS`/`PERIOD` 视需求。

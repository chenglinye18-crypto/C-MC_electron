function [Tprime, p, out] = MC_solve_BTE(cs, mat, opts)
% MC_solve_BTE  Skeleton of a transient deviational MC solver for phonon BTE
%
% USAGE:
%   cs   = setup_case_cube_100nm('Nx',32,'Ny',32,'Nz',32);      % 你的几何/边界
%   si   = mat_silicon_100();                                   % 你的材料/色散
%   opts = mc_default_opts();                                    % 计算参数
%   [Tp, q, out] = MC_solve_BTE(cs, si, opts);
%
% INPUT:
%   geom  struct  几何&网格&边界（通用：立方体/其他几何以后也能复用）
%          .L = [Lx Ly Lz]           尺寸 (m)
%          .Nx, .Ny, .Nz             规则网格计数（便于构建索引/演示）
%          .bc.xmin/xmax/ymin/...    边界: struct('type','isothermal|adiabatic|periodic','T',val,'N',particles)
%          你也可以在 geom 里附带更通用的网格（cell centers/volumes/faces），
%          骨架里留了接口：init_mesh_from_geom_(geom)
%
%   mat   struct  材料/色散（如 mat_silicon_100 的返回）
%          .omega(b,q) .vg(b,q) .energy_meV(b,q) .branch_names, .degeneracy, ...
%
%   opts  struct  运行参数（时间步、粒子数、收敛判据等），可用 mc_default_opts() 生成并修改
%
% OUTPUT:
%   Tprime  (Nc×1)  每个控制体的偏温 T' (K)
%   qfield  struct  各向量通量分量：qfield.qx, qfield.qy, qfield.qz  (Nc×1)  [W/m^2]
%   out     struct  其它信息：时间序列、诊断量等（以后逐步往里加）

% ---------- 0) 参数/随机种子 ----------
if nargin<3 || isempty(opts), opts = mc_default_opts(); end
rng(opts.mc_seed,'twister');

% ---------- 1) 网格与 O(1) 单元定位辅助 ----------
mesh = init_mesh_from_geom_(cs);        % TODO: 填 mesh.centers, .vol, .boxes, .lut 等

% ---------- 2) 构建谱离散（q-网格/模态权重/Cv/τ等） ----------
spec = build_spectral_grid_(mat, opts);   % TODO: 返回 spec.B, spec.q, spec.wq, spec.Cv, spec.vg, spec.tauU, spec.tauI, ...

%%按声子数均匀分配
% ---------- 3) 初始化粒子 ----------
%state = init_state_(mesh, spec, opts);   
%plot_init_diagnostics_(spec, state, mesh, struct('out_dir','fig_init'));

% ---------- 4) 运行MC过程 ----------
%[Tprime, p, out] = MC_time_loop_BTE(mesh, spec, opts, state);

%%按能量均匀分配

state = init_state_Energy(mesh, spec, opts);


[Tprime, p, out] = MC_time_loop_BTE(mesh, spec, opts, state);


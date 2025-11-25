function opts = mc_default_opts()
% 计算参数默认值（你可以在外面拿到这个 struct 后覆盖任意项）
opts.T0      = 350.001;        % 线化点/背景温度 [K]
%opts.T_init = 375;
opts.dt      = 1e-14;        % 初始/基准步
opts.dt_min  = 1e-15;
opts.dt_max  = 1e-11;
opts.t_end   = 5e-10;
opts.fly_mode = 'cell';

opts.use_dynamic_dt = true;
opts.dt_safety_cfl  = 0.5;
opts.p_target       = 0.05;  % 期望每步碰撞概率上限

% 粒子数
opts.mc_face_particles_per_step  = 2e4;
opts.mc_vol_particles_per_step   = 2e4;
opts.mc_scatt_particles_per_step = 2e4;
opts.Neff = 1e3; %计算粒子数 = 实际粒子数/Neff
%如按能量分配
opts.E_eff = 1e-18;

% 统计/稳态
opts.stop_when_steady   = true;
opts.steady_tol_inf     = 1e-2;  % K
opts.steady_tol_l2      = 1e-2;  % 相对
opts.steady_min_steps   = 50;
opts.steady_min_time    = 0.0;
opts.steady_streak_need = 3;

% 热源
opts.enable_volume_Q = 0;

% 随机
opts.mc_seed = uint64(20240511);
opts.use_common_random_numbers = true;
opts.surface_scatter_range = true;

% 谱离散
opts.n_q = 5000;     % 每支的 q 栅格数（Γ->X）
opts.n_w = 1000;
opts.weight_by_Cv_for_Q = true; % 体源谱分配策略示意

% 输出控制
opts.output_every = 10;

% 散射参数
opts.scatter_on = 1;
opts.PP_BL = 1.18e-24;
opts.PP_BTN = 10.5e-13;
opts.PP_BTU = 2.89e-18;
opts.PP_BLTO = 1/3.5e-12;
opts.PB_Tsi = 100e-9;

% U_T关系
opts.T_table_min = 1;         % K
opts.T_table_max = 2000;      % 或更高视材料而定
opts.T_table_n   = 256;       % 网格密度
opts.T_table_log = true;      % 低温对数取点更平滑
opts.invert_Newton_iters = 2; % 牛顿修正步数

% 控制收敛
opts.max_steps = 50000;
opts.conv_tol_inf  = 1e-4;   % K
opts.conv_tol_l2   = 1e-5;   % K
opts.conv_tol_rel  = 1e-6;   % （可选）相对阈值
opts.conv_n_consec = 3;      % 连续满足步数
opts.T_underrelax  = 0.5;    % 欠松弛（有振荡就降一点）

if isempty(gcp('nocreate'))
  parpool('threads');   % R2023a+ 推荐；或者 parpool('local');
end
opts.parallel.use_parfor = true;   % 让 particle_fly_ 选择 parfor 路径
opts.log.on = true;
opts.log.fly_verbose = true;       % 继续保留块/汇总打印

opts.viz.enable    = true;                 % 打开快照
opts.viz.out_dir   = 'viz_snapshots';      % 根目录
opts.viz.run_tag   = datestr(now,'yyyymmdd_HHMMSS');  % 子目录
opts.viz.every_n   = 1;                    % 每步都保存（可改成 5、10 减少量）
opts.viz.colormap  = 'parula';             % 或 'turbo','hot'
% opts.viz.clim    = [300 800];            % 固定色标（不设则自动缩放）

%% ===== options: 使用偏差法=====
opts.deviational = true;           % 打开偏差MC
opts.Tref        = 350;          % 参考温度 T_ref（建议 (Th+Tc)/2）
opts.E_eff       = 1e-19;          % 单个偏差包的“绝对”能量
opts.use_bin_center_w = true;


opts.mode = 'deviational';

%绘图
opts.max_points = 5e4;
opts.save_path='figs/';
end
function si = mat_silicon_100()
% MAT_SILICON_100  Silicon (100) —— 离散化的 q–ω 与 q–v_g 表
% 仅输出离散表 + 查表插值句柄（不暴露解析公式）。
%
% 输入:
%   nq_tab  （可选）q 表格节点数，默认 4097（含端点）
%
% 输出 si 结构:
%   .name, .a0, .qmax
%   .branch_names = {'LA','TA','LO','TO'}
%   .degeneracy   = [1 2 1 2]
%   .B            = 4
%   .q            1×M    （q 网格，含 0 与 qmax）  [1/m]
%   .omega_tab    B×M    （各分支 ω(q)           [rad/s]）
%   .vg_tab       B×M    （各分支 v_g(q)         [m/s]）
%   .omega(b,q)   句柄：对离散表做插值（linear, clamp 到[0,qmax]）
%   .vg(b,q)      句柄：同上
%   .omega_all(q) 返回 B×N 的插值矩阵
%   .vg_all(q)    返回 B×N 的插值矩阵
%   .energy_meV(b,q)  ħω 插值后转 meV

  nq_tab = 5000;

  % —— 基础常量与晶格 —— 
  si.name = 'Silicon (100) — table-driven dispersion';
  si.a0   = 5.431e-10;            % 晶格常数 [m]
  si.qmax = 2*pi/si.a0;           % |q| 的上界（Γ→X）[1/m]

  si.branch_names = {'LA','TA','LO','TO'};
  si.degeneracy   = [1, 2, 1, 2];
  si.B = numel(si.branch_names);

  % —— 文献拟合参数（只用来“离线”生成表；对外不暴露解析式）——
  omega0 = [0.00, 0.00, 9.88, 10.20] * 1e13;   % [rad/s]
  vs     = [9.01, 5.23, 0.00, -2.57] * 1e3;    % [m/s]
  cpar   = [-2.00, -2.26, -1.60,  1.11] * 1e-7;% [m^2/s]

  % —— q 网格（含端点），生成 B×M 的 ω(q) 与 v_g(q) 表 —— 
  q = linspace(0, si.qmax, nq_tab+1);
  B = si.B; M = numel(q);

  omega_tab = zeros(B, M);
  vg_tab    = zeros(B, M);

  for b = 1:B
    % 离线生成表：ω(q)=ω0 + v_s q + c q^2；v_g(q)=dω/dq=v_s+2 c q
    w0 = omega0(b); vs_b = vs(b); c_b = cpar(b);
    omega_tab(b,:) = w0 + vs_b.*q + c_b.*(q.^2);
    omega_tab(b,:) = max(omega_tab(b,:), 0);            % 物理上 ω≥0
    vg_tab(b,:)    = vs_b + 2*c_b.*q;
  end

  % —— 小的数值清理（可选）：避免极小负速度/边界毛刺 —— 
  tiny = 1e-12;
  vg_tab( abs(vg_tab) < tiny ) = 0;

  % —— 对外仅暴露“表 + 查表插值” —— 
  si.q         = q;
  si.omega_tab = omega_tab;
  si.vg_tab    = vg_tab;

  % 插值工具：clamp 到域内，再线性插值；超界时外推=端点常值
  si.omega = @(b,qq) interp_from_table_(qq, q, omega_tab(b,:));
  si.vg    = @(b,qq) interp_from_table_(qq, q, vg_tab(b,:));

  si.omega_all = @(qq) cell2mat( arrayfun(@(b) si.omega(b,qq), 1:B, 'uni',0).' );
  si.vg_all    = @(qq) cell2mat( arrayfun(@(b) si.vg(   b,qq), 1:B, 'uni',0).' );

  % ħ（meV·s）
  hbar_meVs = 6.582119569e-13;
  si.energy_meV = @(b,qq) hbar_meVs .* si.omega(b,qq);
end

% ===== 内部：查表插值（线性 + clamp）=====
function y = interp_from_table_(qq, qtab, ytab)
  % clamp 到 [qmin,qmax]
  qmin = qtab(1); qmax = qtab(end);
  qqc = min(max(qq, qmin), qmax);
  % 线性插值；外推用端点值（更稳）
  y = interp1(qtab, ytab, qqc, 'linear');
  % 缺失（极端数值）回退端点
  if any(isnan(y), 'all')
    y(isnan(y) & qqc<=qmin) = ytab(1);
    y(isnan(y) & qqc>=qmax) = ytab(end);
  end
end

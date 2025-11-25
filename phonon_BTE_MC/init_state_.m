function state = init_state_(mesh, spec, opts)
% INIT_STATE_  论文式初始化：全域随机撒点 + (ω→branch)两步抽样 + 查表 q,|vg|
% 依赖表：spec.q_tbl, spec.w_tbl{b}, spec.vg_tbl{b}
% 依赖频谱：spec.N_w(B×Nw)  (若已含简并度，请直接用；否则见注释)

  if nargin < 3, opts = struct(); end
  if ~isfield(opts,'Neff'),             opts.Neff = 1;        end
  if ~isfield(opts,'use_bin_center_w'), opts.use_bin_center_w = true; end
  if ~isfield(opts,'enhance_factor'),   opts.enhance_factor   = 1;     end
  if ~isfield(opts,'max_particles'),    opts.max_particles    = 2e7;   end

  if isfield(opts,'mc_seed'); rng(opts.mc_seed);
  elseif isfield(opts,'rng_seed'); rng(opts.rng_seed);
  else; rng('shuffle'); end

  % ---- 常量与维度 ----
  hbar = 1.054571817e-34;
  [B,Nw] = size(spec.N_w);
  Nc     = mesh.Nc;
  Vc     = mesh.vol(:);
  Vdom   = sum(Vc);
  Neff   = max(1, round(opts.Neff));
  %g = spec.si.degeneracy(:);  

  have_w_edges = isfield(spec,'w_edges') && numel(spec.w_edges)==(Nw+1);

  % ---- 总数密度 & 能量密度 ----
  Nw_BM          = spec.N_w;                     % B×Nw（数密度×Δω）
  Nw_tot_density = spec.N_w_tot;                 % m^-3
  U_tot_density  = hbar * sum(spec.w_mid .* Nw_BM, 'all');   % J/m^3

  % ---- 总粒子数（随机全域撒点）----
  N_real_total = Nw_tot_density * Vdom;
  Nsp_tot      = round(N_real_total / Neff);
  if Nsp_tot > opts.max_particles
      error('init_state_: 期望粒子数 %g 超过上限 %g；请增大 Neff 或减小体系/温度。',...
            Nsp_tot, opts.max_particles);
  end
  if Nsp_tot==0
      state = empty_state_(); return;
  end

  % ===================== 新：一次抽 (b,m) 的联合 CDF =====================
  % 权重矩阵：W_bm = g_b * N_w(b,m)，再拉平成向量
  W_bm = ones(1,Nw) .* Nw_BM;            % B×Nw
  wvec = reshape(W_bm, [], 1);                  % (B*Nw)×1
  if ~any(wvec>0), error('spec.N_w 全零，无法初始化。'); end
  cdf_bm = cumsum(wvec) / sum(wvec);

  % 工具：从 flat 索引 -> (b,m)
  ind2bm = @(idx) deal( floor((idx-1)/Nw) + 1, 1 + mod(idx-1, Nw) );

  % ===================== 随机把 Nsp_tot 个粒子分到各 cell =====================
  af = opts.enhance_factor; if isscalar(af), af = repmat(af, Nc,1); end
  W_cell   = af(:) .* Vc;                       % 体积×增强因子
  cdf_cell = cumsum(W_cell) / sum(W_cell);
  cid_all  = discretize(rand(Nsp_tot,1), [0; cdf_cell]);    % Nsp_tot×1 ∈[1..Nc]
  Nsp_cell = accumarray(cid_all, 1, [Nc,1], @sum, 0);

  % ===================== 预分配粒子表 =====================
  p(Nsp_tot,1) = struct('id',0,'cell',0, ...
                        'x',0,'y',0,'z',0, ...
                        'b',0,'m',0,'w',0,'q',0, ...
                        'v',zeros(1,3),'vabs',0, ...
                        'E',0,'Wp',Neff, ...
                        'seed',0,'par_id',0,'t_left',0);

  % 是否有全局频率边界（统一分箱）
  have_w_edges = isfield(spec,'w_edges') && numel(spec.w_edges)==(Nw+1);

  % ===================== 主循环：位置→(b,m)→ω→q,|vg|→方向 =====================
  for i = 1:Nsp_tot
      cid = cid_all(i);
      box = mesh.boxes(cid,:);                   % 用该 cell 的盒子，不要用 domain_box
      xmin=box(1); xmax=box(2);
      ymin=box(3); ymax=box(4);
      zmin=box(5); zmax=box(6);

      % 位置：本 cell 内均匀
      x = xmin + (xmax-xmin)*rand();
      y = ymin + (ymax-ymin)*rand();
      z = zmin + (zmax-zmin)*rand();

      % —— 一次抽 (b,m) ——（联合CDF）
      idx_flat = sample_from_cdf_(cdf_bm);
      [b, m]   = ind2bm(idx_flat);

      % bin 内 ω：中心 or 区间均匀
      if opts.use_bin_center_w || ~have_w_edges
          w = spec.w_mid(b,m);                   % 各支一样，这里取 b,m 也OK
      else
          wlo = spec.w_edges(m);  whi = spec.w_edges(m+1);
          w   = wlo + (whi - wlo) * rand();
      end

      % 由表反解 q(w) 并查 |vg(w)|（不使用公式）
      [q, vabs] = q_vabs_from_w_table_(w, spec, b);

      % 各向同性方向（2自由度）
      [ex,ey,ez] = rand_unit_vec_();
      vvec = vabs * [ex,ey,ez];

      % 超粒子能量
      E_sp = Neff * hbar * w;

      % 写入
      p(i).id=i; p(i).cell=cid;
      p(i).x=x;  p(i).y=y;  p(i).z=z;
      p(i).b=b;  p(i).m=m;  p(i).w=w; p(i).q=q;
      p(i).v=vvec; p(i).vabs=vabs;
      p(i).E=E_sp; p(i).Wp=Neff;
      p(i).seed=randi(2^31-1);
      p(i).par_id=i; p(i).t_left=0;
  end

  % ---- 聚合量 ----
  E_current = accumarray([p.cell].', [p.E].', [Nc 1], @sum, 0);
  E_target  = U_tot_density .* Vc;

  state = struct();
  state.p            = p;
  state.Wp           = Neff;
  state.Neff         = Neff;
  state.Nsp_cell     = Nsp_cell;
  state.E_target     = E_target;
  state.E_current    = E_current;
  state.enhance_factor = af;
  state.info = struct('Nsp_tot', Nsp_tot, ...
                      'num_cells', Nc, ...
                      'N_real_total', N_real_total, ...
                      'Nw_tot_density', Nw_tot_density, ...
                      'U_tot_density',  U_tot_density);
end

% ---------- 小工具 ----------
function S = empty_state_()
  S = struct('p',[],'Wp',[],'Neff',[],'Nsp_cell',[],'E_target',[], ...
             'E_current',[],'enhance_factor',[],'info',struct());
end

function idx = sample_from_cdf_(cdf)
  r = rand();
  idx = find(cdf >= r, 1, 'first');
  if isempty(idx), idx = numel(cdf); end
end

function [ex,ey,ez] = rand_unit_vec_()
  u1 = rand(); u2 = rand();
  cz = 2*u1 - 1;
  sz = sqrt(max(0,1 - cz^2));
  phi = 2*pi*u2;
  ex = sz*cos(phi); ey = sz*sin(phi); ez = cz;
end

function [q, vabs] = q_vabs_from_w_table_(w, spec, b)
% q_vabs_from_w_table_  从表反解 q(w) 并取 |vg(w)|
% 输入:
%   w    : 标量或向量 (rad/s)
%   spec : 需含 spec.si.q [Nq×1], spec.si.omega_tab [Nq×B], spec.si.vg_tab [Nq×B]
%   b    : 分支索引
% 输出:
%   q    : 与 w 对应的波矢（与 w 尺寸一致）
%   vabs : |vg|（与 w 尺寸一致）

  qv = spec.si.q(:);                 % Nq×1
  wv = spec.si.omega_tab(:,b);       % Nq×1
  gv = spec.si.vg_tab(:,b);          % Nq×1

  % --- 确保 ω 单调：按 ω 排序，再去重 ---
  [w_sorted, Is] = sort(wv, 'ascend');
  q_sorted = qv(Is);
  v_sorted = gv(Is);

  % 去重（防止重复 ω 点导致 interp1 出错）
  [ws, Iu]  = unique(w_sorted, 'stable');
  qs        = q_sorted(Iu);
  vs        = v_sorted(Iu);

  % --- clamp 到数据范围 ---
  w = w(:);                                    % 向量化
  w_cl = min(max(w, ws(1)), ws(end));

  % --- 插值：q(w) 与 vg(q(w)) ---
  % 若你更偏好直接 vg(ω)，也可：v = interp1(ws, vs, w_cl, 'pchip');
  q    = interp1(ws, qs, w_cl, 'pchip');       % 单调/平滑且不超射
  v    = interp1(qs, vs, q,    'pchip');
  vabs = abs(v);

  % 形状恢复为输入形状
  q    = reshape(q,    size(w_cl));
  vabs = reshape(vabs, size(w_cl));
end

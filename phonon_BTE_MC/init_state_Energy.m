function state = init_state_Energy(mesh, spec, opts)
% INIT_STATE_ENERGY  固定能量超粒子初始化（absolute / deviational）
% opts:
%   .mode  = 'absolute' | 'deviational'   (default 'absolute')
%   .T_init(absolute)                     初始温度
%   .T0, .Tref (deviational)              偏差相对与参考温度
%   .E_eff (default 1e-18)                每个超粒子携带的能量
%   .use_bin_center_w (default true)
%   .enhance_factor (default 1; 可 Nc×1)
%   .max_particles (default 2e8)

  if nargin<3, opts = struct(); end
  if ~isfield(opts,'mode'),               opts.mode = 'absolute'; end
  if ~isfield(opts,'E_eff'),              opts.E_eff = 1e-18; end
  if ~isfield(opts,'use_bin_center_w'),   opts.use_bin_center_w = true; end
  if ~isfield(opts,'enhance_factor'),     opts.enhance_factor = 1; end
  if ~isfield(opts,'max_particles'),      opts.max_particles = 2e8; end

  % ---- 常量 & 频谱表 ----
  kB   = 1.380649e-23;
  hbar = 1.054571817e-34;

  DOSb = spec.DOS_w_b;     % B×Nw  (已含简并度)
  wmid = spec.w_mid;       % B×Nw
  dw   = spec.dw;          % 1×Nw
  [B,Nw] = size(DOSb);

  % ---- 网格体积 ----
  [Nc, Vc] = local_cell_volumes_(mesh);   % Nc×1
  Vdom = sum(Vc);

  % ---- Bose 占据 ----
  bose = @(w,T) 1 ./ max(exp(min(hbar.*w./(kB*T),700)) - 1, realmin);

  % ---- 构造能量密度谱 ----
  switch lower(opts.mode)
    case 'absolute'
      if ~isfield(opts,'T_init'), error('init_state_Energy: 需提供 opts.T_init (absolute)'); end
      nT   = bose(wmid, opts.T_init);                   % B×Nw
      NwT  = DOSb .* nT .* repmat(dw,B,1);              % 每体积个数谱
      Wbm  = hbar .* wmid .* NwT;                       % 绝对能量密度谱 [J/m^3] (>=0)
      Wbm_abs = Wbm;  sgn_bm = ones(B,Nw);
      U_density = sum(Wbm(:));                          % 绝对能量密度
      U_total   = U_density * Vdom;

    case 'deviational'
      if ~isfield(opts,'T0') || ~isfield(opts,'Tref')
        error('init_state_Energy: 需提供 opts.T0 与 opts.Tref (deviational)');
      end
      n0   = bose(wmid, opts.T0);
      nr   = bose(wmid, opts.Tref);
      dn   = n0 - nr;                                   % 可正可负
      Ndev = DOSb .* dn .* repmat(dw,B,1);
      Wbm  = hbar .* wmid .* Ndev;                      % 偏差能量密度谱 (可正可负)
      Wbm_abs = abs(Wbm);  sgn_bm = sign(Wbm); sgn_bm(sgn_bm==0) = 1;
      U_density = sum(Wbm(:));                          % 净偏差能量密度
      U_total   = U_density * Vdom;

    otherwise
      error('opts.mode 仅支持 absolute / deviational');
  end

  % ---- 期望粒子数 ----
  E_eff = opts.E_eff;
  Nexp_tot = (sum(Wbm_abs(:)) * Vdom) / E_eff;
  Nsp_tot  = floor(Nexp_tot) + (rand < (Nexp_tot - floor(Nexp_tot)));

  if Nsp_tot==0
      state = local_empty_state_();
      state.info = struct('mode',opts.mode,'U_total',U_total, ...
                          'Nexp_tot',Nexp_tot,'Nsp_tot',0);
      return;
  end
  if Nsp_tot > opts.max_particles
      error('init_state_Energy: 期望粒子数 %.3g 超上限 %.3g；请增大 E_eff 或缩小系统。', ...
            Nsp_tot, opts.max_particles);
  end

  % ---- (b,m) CDF（按 |Wbm| ）----
  wvec = reshape(Wbm_abs, [], 1);
  cdf_bm = cumsum(wvec) / sum(wvec);
  %ind2bm = @(idx) deal(floor((idx-1)/Nw)+1, 1+mod(idx-1,Nw));
  ind2bm = @(idx) deal(1+mod(idx-1, B), floor((idx-1)/B)+1);

  % ---- 按体积(×增强因子)分配到 cell ----
  af = opts.enhance_factor; if isscalar(af), af = repmat(af,Nc,1); end
  Wcell   = Vc(:) .* af(:);
  cdf_cell= cumsum(Wcell)/sum(Wcell);
  sample_from_cdf = @(cdf) find(cdf >= rand(), 1, 'first');
  cid_all = arrayfun(@(~) sample_from_cdf(cdf_cell), 1:Nsp_tot).';

  % ---- 生成粒子 ----
  have_w_edges = isfield(spec,'w_edges') && numel(spec.w_edges)==(Nw+1);
  p(Nsp_tot,1) = local_blank_particle_();
  for i = 1:Nsp_tot
      cid = cid_all(i);
      [x,y,z] = local_uniform_position_in_cell_(mesh, cid);

      % 选 (b,m)
      idx = local_sample_from_cdf_(cdf_bm);
      [b,m] = ind2bm(idx);

      % ω：中心或 bin 内均匀
      if opts.use_bin_center_w || ~have_w_edges
          w = wmid(b,m);
      else
          w = spec.w_edges(m) + (spec.w_edges(m+1)-spec.w_edges(m))*rand();
      end

      % 从材料表反解 q 与 |vg|
      [q, vabs] = local_q_vabs_from_w_table_(w, spec, b);

      % 方向（体源：整球各向同性）
      dir = local_rand_unit_vec_();
      v   = vabs * dir;

      % 能量符号（偏差）/ 正能量（绝对）
      if strcmpi(opts.mode,'deviational')
          sgn = sgn_bm(b,m);  E = sgn * E_eff;
      else
          sgn = +1;           E = E_eff;
      end

      % 写入
      p(i).id=i; p(i).cell=cid;
      p(i).x=x; p(i).y=y; p(i).z=z;
      p(i).b=b; p(i).m=m; p(i).w=w; p(i).q=q;
      p(i).v=v; p(i).vabs=vabs;
      p(i).E=E; p(i).sgn=sgn;
      p(i).n_ph = E / (hbar*max(w,1e-30));
      p(i).seed=randi(2^31-1); p(i).par_id=i; p(i).t_left=0;
  end

  % ---- 汇总 ----
  Nsp_cell = accumarray([p.cell].', 1, [Nc 1], @sum, 0);
  state = struct();
  state.p = p;
  state.WE = E_eff; state.Wp = E_eff;
  state.Nsp_cell = Nsp_cell;
  state.enhance_factor = af;
  state.info = struct('mode',opts.mode, ...
                      'U_density',U_density, 'U_total',U_total, ...
                      'Nexp_tot',Nexp_tot, 'Nsp_tot',Nsp_tot, ...
                      'Nc',Nc,'Vdom',Vdom);
end

% ----------------- 工具函数 -----------------
function S = local_empty_state_()
  S = struct('p',[],'WE',[],'Wp',[],'Nsp_cell',[], ...
             'E_target',[],'E_current',[],'enhance_factor',[], ...
             'info',struct());
end

function [Nc, Vc] = local_cell_volumes_(mesh)
  if isfield(mesh,'cell_vol') && ~isempty(mesh.cell_vol)
      Vc = mesh.cell_vol(:); Nc = numel(Vc); return;
  end
  if all(isfield(mesh,{'x_edges','y_edges','z_edges','Nx','Ny','Nz'}))
      dx = diff(mesh.x_edges); dy = diff(mesh.y_edges); dz = diff(mesh.z_edges);
      [DX,DY,DZ] = ndgrid(dx,dy,dz);
      Vc = DX(:).*DY(:).*DZ(:); Nc = numel(Vc); return;
  end
  if isfield(mesh,'boxes')
      b = mesh.boxes;
      Vc = (b(:,2)-b(:,1)).*(b(:,4)-b(:,3)).*(b(:,6)-b(:,5));
      Nc = size(b,1); return;
  end
  error('mesh 未提供体积信息。');
end

function [x,y,z] = local_uniform_position_in_cell_(mesh, cid)
  epsl = 1e-12;
  if all(isfield(mesh,{'Nx','Ny','Nz','x_edges','y_edges','z_edges'}))
      [ix,iy,iz] = ind2sub([mesh.Nx,mesh.Ny,mesh.Nz], cid);
      X=mesh.x_edges; Y=mesh.y_edges; Z=mesh.z_edges;
      x = X(ix) + (X(ix+1)-X(ix))*rand();
      y = Y(iy) + (Y(iy+1)-Y(iy))*rand();
      z = Z(iz) + (Z(iz+1)-Z(iz))*rand();
      x = min(max(x,X(ix)+epsl), X(ix+1)-epsl);
      y = min(max(y,Y(iy)+epsl), Y(iy+1)-epsl);
      z = min(max(z,Z(iz)+epsl), Z(iz+1)-epsl);
  elseif isfield(mesh,'boxes')
      b = mesh.boxes(cid,:);
      x = b(1)+(b(2)-b(1))*rand();
      y = b(3)+(b(4)-b(3))*rand();
      z = b(5)+(b(6)-b(5))*rand();
      x = min(max(x,b(1)+epsl), b(2)-epsl);
      y = min(max(y,b(3)+epsl), b(4)-epsl);
      z = min(max(z,b(5)+epsl), b(6)-epsl);
  else
      error('网格信息不足，无法采样位置');
  end
end

function idx = local_sample_from_cdf_(cdf)
  r = rand(); idx = find(cdf>=r,1,'first'); if isempty(idx), idx = numel(cdf); end
end

function u = local_rand_unit_vec_()
  u1 = rand(); u2 = rand(); cz = 2*u1-1;
  sz = sqrt(max(0,1-cz^2)); phi=2*pi*u2; u=[sz*cos(phi), sz*sin(phi), cz];
end

function [q, vabs] = local_q_vabs_from_w_table_(w, spec, b)
% 从材料表反解 q(w) 并取 |vg(q)|
  qv = spec.si.q(:);             % Nq×1
  wv = spec.si.omega_tab(:,b);   % Nq×1
  gv = spec.si.vg_tab(:,b);      % Nq×1

  [w_sorted, Is] = sort(max(wv,0), 'ascend');
  q_sorted = qv(Is); v_sorted = gv(Is);
  [ws, Iu] = unique(w_sorted, 'stable');
  qs = q_sorted(Iu); vs = v_sorted(Iu);

  w_cl = min(max(w, ws(1)), ws(end));
  q    = interp1(ws, qs, w_cl, 'pchip');
  v    = interp1(qs, vs, q,    'pchip');
  vabs = abs(v);
end

function s = local_blank_particle_()
% 统一的粒子字段模板，便于预分配
  s = struct( ...
      'id',0, ...
      'cell',0, ...
      'x',0,'y',0,'z',0, ...
      'b',0,'m',0,'w',0,'q',0, ...
      'v',[0 0 0],'vabs',0, ...
      'E',0,'sgn',+1,'n_ph',0, ...
      'seed',0,'par_id',0,'t_left',0);
end


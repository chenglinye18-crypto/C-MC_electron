function newp = spawn_heat_source(opts, mesh, spec, state, Tprime, LUT, LUTq, src, dt)
% SPAWN_HEAT_SOURCE  体/面热源（谱基准=局部温度 Tloc；查表优先）
% 用法（建议在主循环里先构好 LUT/LUTq，再传进来）：
%   LUT  = build_E_T_lookup_(spec, ...);   % U(T) 与 inv
%   LUTq = build_q_T_lookup_(spec, ...);   % q''(T) 与 inv / inv_dev
%   newp = spawn_heat_source(mesh, spec, state, Tprime, LUT, LUTq, src, dt)
%
% src 需要字段：
%   .type='volume'|'surface';  .qvol 或 .qsurf； .E_eff；
%   体源：src.region 可选 {'cells','box','custom'}（见 region_sampler_volume_）
%   面源：src.face.normal, src.face.bounds（矩形面示例）
%
% 总能量：
%   体源 ΔE = qvol * V_region * dt
%   面源 ΔE = qsurf * A_face   * dt
%
% 反解与谱：
%   体源： U(Tsrc)-U(Tloc) = ΔE/V  → 用 LUT.inv 反解 Tsrc
%   面源： q''(Tsrc)-q''(Tloc) = qsurf → 用 LUTq.inv_dev(Tloc) 反解 Tsrc
%   Wbm ∝ ħω·DOS·[nBE(Tsrc)-nBE(Tloc)]·dw
%
% 发射：
%   先按 |Wbm| 归一化把 |ΔE| 分配到各 (b,m)，再随机四舍五入得到整数粒子数；
%   每个粒子能量 = ±E_eff，符号由 Wbm 的号决定。

  % —— 超粒子能量 E_eff：优先 src.E_eff，否则 opts.E_eff ——
  if ~isfield(src,'E_eff') || isempty(src.E_eff)
      if isfield(opts,'E_eff') && ~isempty(opts.E_eff)
          src.E_eff = opts.E_eff;
      else
          error('spawn_heat_source: 需提供 src.E_eff 或 opts.E_eff');
      end
  end

  % —— 采样用的 Tref（必需）——
  if ~isfield(opts,'Tref') || isempty(opts.Tref)
      if isfield(state,'info') && isfield(state.info,'Tref') && ~isempty(state.info.Tref)
          opts.Tref = state.info.Tref;
      else
          error('spawn_heat_source: 需提供 opts.Tref 或 state.info.Tref 作为采样参考温度');
      end
  end


  switch lower(src.type)
    case 'volume'
      newp = local_spawn_volume_Tloc_(opts, mesh, spec, state, Tprime, LUT, src, dt);
    case 'surface'
      newp = local_spawn_surface_Tloc_(opts,mesh, spec, state, Tprime, LUTq, src, dt);
    otherwise
      error('spawn_heat_source: src.type 仅支持 volume/surface');
  end
end

% ---------------- 体热源：U(Tsrc)-U(Tloc)=ΔE/V ----------------
function p = local_spawn_volume_Tloc_(opts, mesh, spec, state, Tprime, LUT, src, dt)
  E_eff = src.E_eff;

  [cells, Vc, sample_pos_fun, ctr] = region_sampler_volume_(mesh, src);
  V_region = sum(Vc);

  % 局部温度（用 Tprime）
  Tloc = sample_T_at_point_from_Tprime_(mesh, Tprime, ctr);

  % 总能量
  dE_tot = src.qvol * V_region * dt;

  % 反解 Tsrc
  Tsrc = invert_T_from_Udiff_(spec, LUT, Tloc, dE_tot / max(V_region,eps));

  % 局部差分谱
  Wbm = build_local_diff_spectrum_(spec, Tsrc, Tloc);
  if all(Wbm(:)==0), p = struct([]); return; end

  % 发射
  next_id_base = get_next_id_(state);
  p = emit_particles_from_Wbm_refsample_( ...
        mesh, spec, cells, Vc, sample_pos_fun, ...
        Wbm, dE_tot, src.E_eff, [], opts.Tref, next_id_base);
end

% ---------------- 面热源：q''(Tsrc)-q''(Tloc)=qsurf ----------------
function p = local_spawn_surface_Tloc_(opts, mesh, spec, state, Tprime, LUTq, src, dt)
  E_eff = src.E_eff;

  [A_face, sample_on_face, n_hat, ctr] = surface_sampler_(mesh, src.face);

  % 局部温度（用 Tprime）
  Tloc = sample_T_at_point_from_Tprime_(mesh, Tprime, ctr);

  dE_tot = src.qsurf * A_face * dt;

  % 反解 Tsrc
  Tsrc = invert_T_from_qdiff_(spec, LUTq, Tloc, src.qsurf);

  % 局部差分谱
  Wbm = build_local_diff_spectrum_(spec, Tsrc, Tloc);
  if all(Wbm(:)==0), p = struct([]); return; end

  % 发射（半空间余弦）
  next_id_base = get_next_id_(state);
  p = emit_particles_from_Wbm_refsample_( ...
        mesh, spec, [], [], sample_on_face, ...
        Wbm, dE_tot, src.E_eff, n_hat, opts.Tref, next_id_base);

end

% ---------------- 局部差分谱 Wbm ----------------
function Wbm = build_local_diff_spectrum_(spec, Tsrc, Tloc)
  kB=1.380649e-23; hbar=1.054571817e-34;
  bose = @(w,T) 1./max(exp(min(hbar.*w./(kB*T),700))-1, realmin);

  w   = max(spec.w_mid,0);
  DOS = max(spec.DOS_w_b,0);
  if isvector(spec.dw), dw = repmat(reshape(spec.dw,1,[]), size(w,1),1);
  else,                 dw = spec.dw; end

  n_src = bose(w, Tsrc);
  n_loc = bose(w, Tloc);
  Wbm   = hbar .* w .* DOS .* (n_src - n_loc) .* dw;
  Wbm(DOS<=0) = 0;
end

% ---------------- 反解工具（查表优先） ----------------
function Tsrc = invert_T_from_Udiff_(spec, LUT, Tloc, dU)
  if ~isempty(LUT) && isfield(LUT,'inv') && isfield(LUT,'T') && isfield(LUT,'U')
      Uloc = interp1(LUT.T, LUT.U, Tloc, 'pchip','extrap');
      Tsrc = LUT.inv(Uloc + dU);
  else
      Tsrc = invert_monotone(@(T) U_density_equil(spec,T), U_density_equil(spec,Tloc)+dU, 1, 5000);
  end
end

function Tsrc = invert_T_from_qdiff_(spec, LUTq, Tloc, qsurf)
  if ~isempty(LUTq) && isfield(LUTq,'inv_dev')
      inv_fun = LUTq.inv_dev(Tloc);   % 先得到 q→T 的函数句柄
      Tsrc    = inv_fun(qsurf);       % 再用它去解 Tsrc
  else
      vg     = ensure_vg_table(spec);
      target = qflux_equil(spec, Tloc, vg) + qsurf;
      Tsrc   = invert_monotone(@(T) qflux_equil(spec,T,vg), target, 1, 5000);
  end
end


% ---------------- 发射器（按谱分配 |ΔE|，号由谱决定） ----------------
function p = emit_particles_from_Wbm_refsample_(mesh, spec, cells, Vc, sample_pos_fun, ...
                                                Wbm, dE_tot, E_eff, n_hat, Tref, next_id_base)
  % 目标谱(定号)：Wbm ∝ ħω·DOS·(n(Tsrc)-n(Tloc))·dw   [可正可负]
  % 采样谱(定分布)：Wref ∝ ħω·DOS·n(Tref)·dw           [全正，用于抽(b,m,ω)]
  [B,Nw] = size(Wbm);
  Wabs = abs(Wbm); S = sum(Wabs(:));
  if S==0 || ~isfinite(S), p = struct([]); return; end

  % 参考谱（Tref）
  kB=1.380649e-23; hbar=1.054571817e-34;
  w   = max(spec.w_mid,0);
  DOS = max(spec.DOS_w_b,0);
  if isvector(spec.dw), dw = repmat(reshape(spec.dw,1,[]), B,1); else, dw = spec.dw; end
  xref = (hbar.*w)/(kB*max(Tref,1e-12));
  nref = 1 ./ max(exp(min(xref,700))-1, realmin);
  Wref = hbar.*w .* DOS .* nref .* dw;
  if all(Wref(:)==0), p = struct([]); return; end
  cdf_ref = cumsum(Wref(:)) / sum(Wref(:));

  % 仅允许抽到目标谱非零的格
  nonzero_mask = (Wabs>0);
  if ~any(nonzero_mask(:)), p = struct([]); return; end

  % 本步总粒子数
  Nexp = abs(dE_tot)/E_eff;
  Nsp  = floor(Nexp) + (rand < (Nexp - floor(Nexp)));
  if Nsp==0, p = struct([]); return; end

  if ~isempty(cells)
      cdf_cell = cumsum(Vc(:))/sum(Vc(:));
  end

  p(Nsp,1) = local_blank_particle_();
  has_edges = isfield(spec,'w_edges') && numel(spec.w_edges)==(Nw+1);

  % 从参考CDF抽一个非零目标格
  pick_bm = @() local_pick_bm_nonzero_(cdf_ref, B, Nw, nonzero_mask);

  for i = 1:Nsp
      [bb,mm] = pick_bm();
      sgn = sign(Wbm(bb,mm)); if sgn==0, sgn=+1; end

  if ~isempty(cells)
      cid = cells(find(cdf_cell>=rand(),1,'first'));
      [x,y,z] = local_uniform_position_in_cell_(mesh, cid);
  else
      [x,y,z] = sample_pos_fun();
      % 反定位 cell
      cid = locate_cell_from_point_(mesh, [x,y,z]);
      if cid==0
          % 边界上的点，沿法向（若有）或用夹取法微移一丝再定位
          if ~isempty(n_hat)
              epsl = 1e-12;
              nn = n_hat(:).'/max(norm(n_hat),1e-30);
              x = x + epsl*nn(1);  y = y + epsl*nn(2);  z = z + epsl*nn(3);
              cid = locate_cell_from_point_(mesh, [x,y,z]);
          end
          if cid==0
              [x,y,z,cid] = local_nudge_point_inside_(mesh, [x,y,z]);
          end
      end
  end



      if has_edges
          w_i = spec.w_edges(mm) + (spec.w_edges(mm+1)-spec.w_edges(mm))*rand();
      else
          w_i = w(bb,mm);
      end
      [q_i, vabs_i] = local_q_vabs_from_w_table_(w_i, spec, bb);

      if isempty(n_hat)
          dir = local_rand_unit_vec_();
      else
          nh = n_hat(:).'/max(norm(n_hat),1e-30);
          mu = sqrt(rand()); phi=2*pi*rand();
          u_local = [sqrt(1-mu^2)*cos(phi), sqrt(1-mu^2)*sin(phi), mu];
          R = local_rot_from_z_to_n_(nh); dir = (R*u_local.').';
      end
      v_i = vabs_i * dir;

      pid = next_id_base + i;
      P = local_blank_particle_();
      P.id=pid; P.par_id=pid; P.cell=cid; P.x=x; P.y=y; P.z=z;
      P.b=bb; P.m=mm; P.w=w_i; P.q=q_i; P.v=v_i; P.vabs=vabs_i;
      P.E = sgn*E_eff; P.sgn=sgn;
      P.n_ph = (sgn*E_eff)/(hbar*max(w_i,1e-30));
      p(i)=P;
  end
end

function [bb,mm] = local_pick_bm_nonzero_(cdf_ref, B, Nw, mask)
  for tries=1:20
      r = rand(); idx = find(cdf_ref>=r,1,'first'); if isempty(idx), idx=numel(cdf_ref); end
      bb = 1 + mod(idx-1, B);
      mm = floor((idx-1)/B) + 1;
      if mask(bb,mm), return; end
  end
  [bb,mm] = find(mask,1,'first');
end

% ---------------- 区域/面采样（含中心点，用于取 Tloc） ----------------
function [cells, Vc, sample_pos_fun, ctr] = region_sampler_volume_(mesh, src)
  if isfield(src,'region') && isstruct(src.region), R=src.region; else, R=struct('type','cells','id',[]); end
  switch lower(R.type)
    case 'cells'
      [Nc,Vc_all] = local_cell_volumes_(mesh);
      if ~isfield(R,'id') || isempty(R.id), cells=(1:Nc).'; else, cells=R.id(:); end
      Vc = Vc_all(cells);
      cdf = cumsum(Vc)/sum(Vc);
      sample_pos_fun = @() local_uniform_position_in_cell_(mesh, cells(find(cdf>=rand(),1,'first')));
      ctr = estimate_cells_center_(mesh, cells);
    case 'box'
      b = R.bounds; Vc = (b(2)-b(1))*(b(4)-b(3))*(b(6)-b(5)); cells=[];
      sample_pos_fun = @() deal(b(1)+(b(2)-b(1))*rand(), b(3)+(b(4)-b(3))*rand(), b(5)+(b(6)-b(5))*rand());
      ctr = [(b(1)+b(2))/2, (b(3)+b(4))/2, (b(5)+b(6))/2];
    case 'custom'
      cells=[]; Vc = R.measure; sample_pos_fun = @() R.sample_pos();
      if isfield(R,'center'), ctr = R.center; else, ctr = sample_pos_fun(); end
    otherwise
      error('region_sampler_volume_: 未知 region.type');
  end
end

function ctr = estimate_cells_center_(mesh, cells)
  if all(isfield(mesh,{'x_edges','y_edges','z_edges','Nx','Ny','Nz'}))
      [ix,iy,iz] = ind2sub([mesh.Nx,mesh.Ny,mesh.Nz], round(mean(cells)));
      ctr = [ mean(mesh.x_edges([ix ix+1])), ...
              mean(mesh.y_edges([iy iy+1])), ...
              mean(mesh.z_edges([iz iz+1])) ];
  elseif isfield(mesh,'boxes')
      b = mesh.boxes(cells,:); c1 = mean(b(:,[1 3 5]),1); c2 = mean(b(:,[2 4 6]),1);
      ctr = 0.5*(c1+c2);
  else
      ctr = [0 0 0];
  end
end

function [A_face, sample_on_face, n_hat, ctr] = surface_sampler_(mesh, face)
  n_hat = face.normal(:).'/norm(face.normal);
  b = face.bounds;
  if abs(n_hat(1))==1       % yz 面 x=const
      A_face = (b(4)-b(3))*(b(2)-b(1));
      sample_on_face = @() deal(b(5), b(3)+(b(4)-b(3))*rand(), b(1)+(b(2)-b(1))*rand());
      ctr = [b(5), (b(3)+b(4))/2, (b(1)+b(2))/2];
  elseif abs(n_hat(2))==1   % xz 面 y=const
      A_face = (b(2)-b(1))*(b(6)-b(5));
      sample_on_face = @() deal(b(1)+(b(2)-b(1))*rand(), b(5), b(3)+(b(4)-b(3))*rand());
      ctr = [(b(1)+b(2))/2, b(5), (b(3)+b(4))/2];
  else                      % xy 面 z=const
      A_face = (b(2)-b(1))*(b(4)-b(3));
      sample_on_face = @() deal(b(1)+(b(2)-b(1))*rand(), b(3)+(b(4)-b(3))*rand(), b(5));
      ctr = [(b(1)+b(2))/2, (b(3)+b(4))/2, b(5)];
  end
end

% ---------------- 从 Tprime 取局部温度 ----------------
function Tloc = sample_T_at_point_from_Tprime_(mesh, Tprime, pt)
  cid = locate_cell_from_point_(mesh, pt);
  if cid>0 && cid<=numel(Tprime), Tloc = Tprime(cid);
  else, Tloc = mean(Tprime); % 出界/未命中时取均值兜底
  end
end

function cid = locate_cell_from_point_(mesh, pt)
  x=pt(1); y=pt(2); z=pt(3);
  if all(isfield(mesh,{'x_edges','y_edges','z_edges','Nx','Ny','Nz'}))
      ix = find(mesh.x_edges<=x,1,'last');
      iy = find(mesh.y_edges<=y,1,'last');
      iz = find(mesh.z_edges<=z,1,'last');
      if isempty(ix)||isempty(iy)||isempty(iz)||ix>=numel(mesh.x_edges)||iy>=numel(mesh.y_edges)||iz>=numel(mesh.z_edges)
          cid = 0; return;
      end
      cid = sub2ind([mesh.Nx,mesh.Ny,mesh.Nz], ix, iy, iz);
  elseif isfield(mesh,'boxes')
      cid=0; bx=mesh.boxes;
      for i=1:size(bx,1)
          if x>=bx(i,1)&&x<bx(i,2)&&y>=bx(i,3)&&y<bx(i,4)&&z>=bx(i,5)&&z<bx(i,6), cid=i; break; end
      end
  else
      cid = 0;
  end
end

% ---------------- 物性工具（若没 LUT 则现算） ----------------
function U = U_density_equil(spec, T)
  kB=1.380649e-23; hbar=1.054571817e-34;
  w=spec.w_mid; DOS=spec.DOS_w_b; 
  if isvector(spec.dw), dw=repmat(reshape(spec.dw,1,[]),size(w,1),1); else, dw=spec.dw; end
  n = 1 ./ max(exp(min(hbar.*w./(kB*T),700))-1, realmin);
  U = sum(sum(hbar.*w.*DOS.*n.*dw));
end

function q = qflux_equil(spec, T, vg)
  kB=1.380649e-23; hbar=1.054571817e-34;
  w=spec.w_mid; DOS=spec.DOS_w_b;
  if isvector(spec.dw), dw=repmat(reshape(spec.dw,1,[]),size(w,1),1); else, dw=spec.dw; end
  if nargin<3||isempty(vg), vg=ensure_vg_table(spec); end
  n = 1 ./ max(exp(min(hbar.*w./(kB*T),700))-1, realmin);
  q = 0.25*sum(sum(DOS.*vg.*(hbar.*w).*n.*dw));
end

function vg = ensure_vg_table(spec)
  [B,Nw]=size(spec.w_mid); vg=zeros(B,Nw);
  for b=1:B, for m=1:Nw
      w=spec.w_mid(b,m); [~,v]=local_q_vabs_from_w_table_(w,spec,b); vg(b,m)=v;
  end, end
end

function Tsol = invert_monotone(fun, target, Tmin, Tmax, tol)
  if nargin<5, tol=1e-6; end
  fmin=fun(Tmin)-target; fmax=fun(Tmax)-target;
  if fmin>0, error('invert_monotone: lower bound too high'); end
  if fmax<0, error('invert_monotone: upper bound too low'); end
  for k=1:100
      Tmid=0.5*(Tmin+Tmax); fmid=fun(Tmid)-target;
      if abs(fmid)<=max(1e-12, tol*max(1,abs(target))), Tsol=Tmid; return; end
      if fmid>0, Tmax=Tmid; else, Tmin=Tmid; end
  end
  Tsol=0.5*(Tmin+Tmax);
end

function nid = get_next_id_(state)
  if isempty(state.p), nid = 0;
  else
      ids = [state.p.id]; if isempty(ids), nid = 0; else, nid = max(ids); end
  end
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

function R = local_rot_from_z_to_n_(n)
% 把 z 轴旋到单位向量 n 的 3×3 旋转矩阵（Rodrigues 公式）
  ez = [0;0;1];
  n  = n(:);
  if norm(n) < 1e-30
      R = eye(3); return;
  end
  n = n / norm(n);
  v = cross(ez, n);
  s = norm(v);
  c = dot(ez, n);
  if s < 1e-12
      % 平行或反平行
      if c > 0
          R = eye(3);
      else
          R = diag([1 -1 -1]);  % 180°绕 x 轴（任意与 z 垂直轴都可）
      end
      return;
  end
  vx = [   0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2)  v(1)   0  ];
  R = eye(3) + vx + vx*vx*((1-c)/(s^2));
end

function [x2,y2,z2,cid] = local_nudge_point_inside_(mesh, pt)
% 把点夹到域内一点点，再做 cell 定位；用于边界/数值误差的兜底
  epsl = 1e-12;
  if all(isfield(mesh,{'x_edges','y_edges','z_edges'}))
      x2 = min(max(pt(1), mesh.x_edges(1)+epsl), mesh.x_edges(end)-epsl);
      y2 = min(max(pt(2), mesh.y_edges(1)+epsl), mesh.y_edges(end)-epsl);
      z2 = min(max(pt(3), mesh.z_edges(1)+epsl), mesh.z_edges(end)-epsl);
  else
      % 没有 edges 信息就微移一丝（不推荐，但保底）
      x2 = pt(1) + epsl; y2 = pt(2) + epsl; z2 = pt(3) + epsl;
  end
  cid = locate_cell_from_point_(mesh, [x2,y2,z2]);
end


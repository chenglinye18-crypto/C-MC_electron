function state = seed_ghost_particles_(state, mesh_g, ghost, spec, opts, bc)
% 在每个 ghost 面对应的 ghost cell 内，按 Tb 的偏差分布采样粒子并追加到 state.p
% 要求：opts.mode='deviational' 且提供 opts.Tref；每个 isothermal 面在 bc.face.T 给出 Tb

  if ~isfield(opts,'mode') || ~strcmpi(opts.mode,'deviational')
    error('seed_ghost_particles_: 请在偏差法下使用（opts.mode=''deviational''）');
  end
  if ~isfield(opts,'Tref'), error('seed_ghost_particles_: 缺少 opts.Tref'); end
  if ~isfield(opts,'E_eff'), opts.E_eff = 1e-18; end
  if ~isfield(opts,'use_bin_center_w'), opts.use_bin_center_w = true; end

  % cell 体积
  Vc = cell_volumes_from_edges_(mesh_g);

  % 频谱准备
  kB=1.380649e-23; hbar=1.054571817e-34;
  DOS=spec.DOS_w_b; wmid=spec.w_mid; vg = abs(spec.vg_w);
  dw = spec.dw; [B,Nw]=size(DOS); if isvector(dw), dw=repmat(reshape(dw,1,[]),B,1); end

  % 参考占据
  xref = (hbar.*wmid)/(kB*max(opts.Tref,1e-12));
  nref = 1./max(exp(min(xref,700))-1, realmin);

  faces = fieldnames(ghost.by_face);
  for f=1:numel(faces)
    tag = faces{f};
    if ~isfield(bc,tag) || ~isfield(bc.(tag),'T'), continue; end
    Tb = bc.(tag).T; if ~(isfinite(Tb) && Tb>0), continue; end

    xTb = (hbar.*wmid)/(kB*Tb);
    nTb = 1./max(exp(min(xTb,700))-1, realmin);

    dN  = nTb - nref;                     % 频谱偏差
    Wbm = (hbar.*wmid).*dN.*dw.*DOS;      % J/m^3（可正可负）
    Wabs= abs(Wbm); sgn = sign(Wbm); sgn(sgn==0)=1;

    % 该 ghost 面的所有 ghost cell IDs
    cids = ghost.by_face.(tag)(:);
    Vsum = sum(Vc(cids));

    % 期望粒子数
    E_eff=opts.E_eff;
    Nexp = (sum(Wabs(:))*Vsum)/E_eff;
    Nadd = floor(Nexp) + (rand < (Nexp - floor(Nexp)));
    if Nadd==0, continue; end

    % (b,m) CDF
    wvec = Wabs(:); cdf_bm = cumsum(wvec)/sum(wvec);
    %ind2bm = @(idx) deal(floor((idx-1)/Nw)+1, 1+mod(idx-1,Nw));
    ind2bm = @(idx) deal(1+mod(idx-1, B), floor((idx-1)/B)+1);

    % 在这组 ghost cells 里按体积均匀选 cell
    cdf_cell = cumsum(Vc(cids))/sum(Vc(cids));
    pick_cell = @(r) cids(find(cdf_cell>=r,1,'first'));

    % 生成
    newP(Nadd,1) = blank_particle_();
    if ~isempty(state.p)
        next_id = state.p(end).id;
    else
        next_id = 1;
    end
    have_w_edges = isfield(spec,'w_edges') && numel(spec.w_edges)==(Nw+1);

    for i=1:Nadd
      cid = pick_cell(rand());
      [x,y,z] = uniform_position_in_cell_(mesh_g, cid);

      idx = find(cdf_bm>=rand(),1,'first'); if isempty(idx), idx=numel(cdf_bm); end
      [b,m] = ind2bm(idx);
      if opts.use_bin_center_w || ~have_w_edges
        w = wmid(b,m);
      else
        w = spec.w_edges(m) + (spec.w_edges(m+1)-spec.w_edges(m))*rand();
      end

      [q, vabs] = q_vabs_from_w_table_(w, spec, b);
      dir = rand_unit_vec_(); v = vabs*dir;
      E   = sign( sgn(b,m) ) * E_eff;

      pid = next_id+i;
      P = blank_particle_();
      P.id=pid; P.cell=cid; P.x=x; P.y=y; P.z=z;
      P.b=b; P.m=m; P.w=w; P.q=q; P.v=v; P.vabs=vabs;
      P.E=E; P.sgn=sign(E); P.n_ph = E/(hbar*w);
      newP(i)=P;
    end

    if isempty(state.p), state.p=newP; else, state.p=[state.p; newP]; end
  end
end

function Vc = cell_volumes_from_edges_(mesh)
  dx = diff(mesh.x_edges(:)); dy = diff(mesh.y_edges(:)); dz = diff(mesh.z_edges(:));
  [DX,DY,DZ] = ndgrid(dx,dy,dz);
  Vc = (DX(:).*DY(:).*DZ(:));
end

function [x,y,z] = uniform_position_in_cell_(mesh, cid)
  [ix,iy,iz] = ind2sub([mesh.Nx,mesh.Ny,mesh.Nz], cid);
  X=mesh.x_edges; Y=mesh.y_edges; Z=mesh.z_edges;
  epsl=1e-12;
  x = X(ix) + (X(ix+1)-X(ix))*rand();
  y = Y(iy) + (Y(iy+1)-Y(iy))*rand();
  z = Z(iz) + (Z(iz+1)-Z(iz))*rand();
  x = min(max(x,X(ix)+epsl), X(ix+1)-epsl);
  y = min(max(y,Y(iy)+epsl), Y(iy+1)-epsl);
  z = min(max(z,Z(iz)+epsl), Z(iz+1)-epsl);
end

function v = getfield_or(s, name, default_v)
  if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name);
  else, v = default_v; end
end

function s = blank_particle_()
  s = struct('id',0,'cell',0,'x',0,'y',0,'z',0,'b',0,'m',0,'w',0,'q',0, ...
             'v',[0 0 0],'vabs',0,'E',0,'sgn',+1,'n_ph',0,'seed',0,'par_id',0,'t_left',0);
end

function [q, vabs] = q_vabs_from_w_table_(w, spec, b)
  qv = spec.si.q(:); wv = spec.si.omega_tab(:,b); gv = spec.si.vg_tab(:,b);
  [w_sorted, Is] = sort(wv, 'ascend');
  q_sorted = qv(Is); v_sorted = gv(Is);
  [ws, Iu]  = unique(w_sorted, 'stable');
  qs = q_sorted(Iu); vs = v_sorted(Iu);
  w_cl = min(max(w, ws(1)), ws(end));
  q    = interp1(ws, qs, w_cl, 'pchip');
  v    = interp1(qs, vs, q,    'pchip');
  vabs = abs(v);
end

function u = rand_unit_vec_()
  u1=rand(); u2=rand(); cz=2*u1-1; sz=sqrt(max(0,1-cz^2)); phi=2*pi*u2;
  u=[sz*cos(phi), sz*sin(phi), cz];
end

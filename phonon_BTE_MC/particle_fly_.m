function [state, absorb] = particle_fly_(state, mesh, dt, opts, spec)
% 纯推进：不做接口能量记账（ghost↔core 的能量统计改为飞行后对比实现）
% mesh.bc: periodic/adiabatic/isothermal/absorb；最外层 absorb 仅移除粒子，不记能量。
  if nargin<4 || isempty(opts), opts = struct(); end
  if ~isfield(opts,'fly_mode'), opts.fly_mode = 'cell'; end
  must = {'Nx','Ny','Nz','x_edges','y_edges','z_edges'};
  if ~all(isfield(mesh,must)), error('particle_fly_: 仅规则网格'); end

  absorb = struct();  % 兼容返回，占位

  switch lower(opts.fly_mode)
    case 'cell',   state = fly_cell_faces_pure_(state, mesh, dt, opts);
    case 'domain', state = fly_domain_only_pure_(state, mesh, dt, opts);
    otherwise, error('particle_fly_: 未知 fly_mode');
  end
end

% ----------------- cell 面推进（纯推进） -----------------
function state = fly_cell_faces_pure_(state, mesh, dt, opts)
  p = state.p; Np=numel(p); if Np==0, return; end

  Nx=mesh.Nx; Ny=mesh.Ny; Nz=mesh.Nz;
  X=mesh.x_edges(:); Y=mesh.y_edges(:); Z=mesh.z_edges(:);
  Xmin=X(1); Xmax=X(end); Ymin=Y(1); Ymax=Y(end); Zmin=Z(1); Zmax=Z(end);
  epsl=1e-11;

  CHUNK = getfield_or(opts,'log',struct()); CHUNK = getfield_or(CHUNK,'fly_chunk',2e5);
  if ~isscalar(CHUNK) || CHUNK<=0, CHUNK = Np; end
  nBlocks = ceil(Np/CHUNK);

  for bId=1:nBlocks
    i1=(bId-1)*CHUNK+1; i2=min(bId*CHUNK,Np); id=(i1:i2).'; nc=numel(id);

    x=[p(id).x].'; y=[p(id).y].'; z=[p(id).z].';
    v=vertcat(p(id).v); vx=v(:,1); vy=v(:,2); vz=v(:,3);
    cid=[p(id).cell].'; vabs=[p(id).vabs].';
    alive=(cid>0) & (vabs>0); t_rem=dt*ones(nc,1);

    [ix,iy,iz] = ind2sub([Nx,Ny,Nz], cid);
    x = min(max(x, X(ix)+epsl), X(ix+1)-epsl);
    y = min(max(y, Y(iy)+epsl), Y(iy+1)-epsl);
    z = min(max(z, Z(iz)+epsl), Z(iz+1)-epsl);

    rel_tol = @(t) max(1e-15*t, 1e-18);

    while true
      act = find(alive & (t_rem>0)); if isempty(act), break; end
      xa=x(act); ya=y(act); za=z(act);
      vxa=vx(act); vya=vy(act); vza=vz(act);
      ix_a=ix(act); iy_a=iy(act); iz_a=iz(act);

      INF=inf(numel(act),1);
      tx=INF; ty=INF; tz=INF;
      xL=X(ix_a); xR=X(ix_a+1);
      yB=Y(iy_a); yT=Y(iy_a+1);
      zD=Z(iz_a); zU=Z(iz_a+1);
      pos=vxa>0; if any(pos), tx(pos)=(xR(pos)-xa(pos))./vxa(pos); end
      neg=vxa<0; if any(neg), tx(neg)=(xL(neg)-xa(neg))./vxa(neg); end
      pos=vya>0; if any(pos), ty(pos)=(yT(pos)-ya(pos))./vya(pos); end
      neg=vya<0; if any(neg), ty(neg)=(yB(neg)-ya(neg))./vya(neg); end
      pos=vza>0; if any(pos), tz(pos)=(zU(pos)-za(pos))./vza(pos); end
      neg=vza<0; if any(neg), tz(neg)=(zD(neg)-za(neg))./vza(neg); end
      tx=max(tx,0); ty=max(ty,0); tz=max(tz,0);

      tcell=tx; axis_ix=ones(numel(act),1,'uint8');
      m=ty<tcell; tcell(m)=ty(m); axis_ix(m)=2;
      m=tz<tcell; tcell(m)=tz(m); axis_ix(m)=3;

      tf=min(tcell, t_rem(act));
      xa=xa+vxa.*tf; ya=ya+vya.*tf; za=za+vza.*tf;
      t_rem(act)=t_rem(act)-tf;

      hit=isfinite(tcell) & (abs(tf-tcell) <= rel_tol(tcell));
      if any(hit)
        ids=find(hit).';
        for jj=ids
          k=act(jj); ax=axis_ix(jj);
          switch ax
            case 1 % X
              if vxa(jj)>0
                if ix(k)>=Nx
                  bc = lower(getfield_or(mesh.bc,'x_max',struct('type','periodic')).type);
                  if strcmp(bc,'periodic')
                    ix(k)=1; xa(jj)=X(1)+epsl;
                  elseif strcmp(bc,'adiabatic')
                    vx(k)=-vx(k); xa(jj)=X(end)-epsl;
                  else
                    alive(k)=false; % absorb 外层：直接丢弃
                  end
                else
                  ix(k)=ix(k)+1; xa(jj)=X(ix(k))+epsl;
                end
              else
                if ix(k)<=1
                  bc = lower(getfield_or(mesh.bc,'x_min',struct('type','periodic')).type);
                  if strcmp(bc,'periodic')
                    ix(k)=Nx; xa(jj)=X(end)-epsl;
                  elseif strcmp(bc,'adiabatic')
                    vx(k)=-vx(k); xa(jj)=X(1)+epsl;
                  else
                    alive(k)=false;
                  end
                else
                  ix(k)=ix(k)-1; xa(jj)=X(ix(k)+1)-epsl;
                end
              end
            case 2 % Y
              if vya(jj)>0
                if iy(k)>=Ny
                  bc = lower(getfield_or(mesh.bc,'y_max',struct('type','periodic')).type);
                  if strcmp(bc,'periodic')
                    iy(k)=1; ya(jj)=Y(1)+epsl;
                  elseif strcmp(bc,'adiabatic')
                    vy(k)=-vy(k); ya(jj)=Y(end)-epsl;
                  else
                    alive(k)=false;
                  end
                else
                  iy(k)=iy(k)+1; ya(jj)=Y(iy(k))+epsl;
                end
              else
                if iy(k)<=1
                  bc = lower(getfield_or(mesh.bc,'y_min',struct('type','periodic')).type);
                  if strcmp(bc,'periodic')
                    iy(k)=Ny; ya(jj)=Y(end)-epsl;
                  elseif strcmp(bc,'adiabatic')
                    vy(k)=-vy(k); ya(jj)=Y(1)+epsl;
                  else
                    alive(k)=false;
                  end
                else
                  iy(k)=iy(k)-1; ya(jj)=Y(iy(k)+1)-epsl;
                end
              end
            otherwise % Z
              if vza(jj)>0
                if iz(k)>=Nz
                  bc = lower(getfield_or(mesh.bc,'z_max',struct('type','periodic')).type);
                  if strcmp(bc,'periodic')
                    iz(k)=1; za(jj)=Z(1)+epsl;
                  elseif strcmp(bc,'adiabatic')
                    vz(k)=-vz(k); za(jj)=Z(end)-epsl;
                  else
                    alive(k)=false;
                  end
                else
                  iz(k)=iz(k)+1; za(jj)=Z(iz(k))+epsl;
                end
              else
                if iz(k)<=1
                  bc = lower(getfield_or(mesh.bc,'z_min',struct('type','periodic')).type);
                  if strcmp(bc,'periodic')
                    iz(k)=Nz; za(jj)=Z(end)-epsl;
                  elseif strcmp(bc,'adiabatic')
                    vz(k)=-vz(k); za(jj)=Z(1)+epsl;
                  else
                    alive(k)=false;
                  end
                else
                  iz(k)=iz(k)-1; za(jj)=Z(iz(k)+1)-epsl;
                end
              end
          end
        end
      end

      x(act)=xa; y(act)=ya; z(act)=za;
      alive_now=find(alive);
      x(alive_now)=min(max(x(alive_now), X(ix(alive_now))+epsl), X(ix(alive_now)+1)-epsl);
      y(alive_now)=min(max(y(alive_now), Y(iy(alive_now))+epsl), Y(iy(alive_now)+1)-epsl);
      z(alive_now)=min(max(z(alive_now), Z(iz(alive_now))+epsl), Z(iz(alive_now)+1)-epsl);
    end

    cid_new = int32(sub2ind([Nx,Ny,Nz], ix,iy,iz));
    for k=1:nc
      ii=id(k);
      if alive(k)
        state.p(ii).x=x(k); state.p(ii).y=y(k); state.p(ii).z=z(k);
        state.p(ii).v=[vx(k) vy(k) vz(k)];
        state.p(ii).cell=cid_new(k); state.p(ii).t_left=0;
      else
        state.p(ii).cell = int32(-1); % 标记已死
      end
    end
  end

  alive_glob = arrayfun(@(pp) pp.cell>0, state.p);
  state.p = state.p(alive_glob);
end

% ----------------- 域边界推进（纯推进） -----------------
function state = fly_domain_only_pure_(state, mesh, dt, opts)
  % 简化：与上相同思路，只在最外层做 periodic/adiabatic/absorb，不再赘述
  % 为保持篇幅精简，建议统一用 cell 模式；若你需要，我再给完整 domain-only 版
  state = fly_cell_faces_pure_(state, mesh, dt, opts);
end

% ----------------- 小工具 -----------------
function v = getfield_or(s, name, default_v)
  if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), v=s.(name); else, v=default_v; end
end

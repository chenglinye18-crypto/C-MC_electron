function iface = tally_iface_by_cid_(p_before, p_after, mesh_g, opts)
% 依据“飞行前/后”的 cell 线性下标，按最简单的一层 ghost 规则统计 ghost↔core 接口能量。
%
% 约定（与打印一致）：
%   E_in  : ghost -> core  的总能量（带符号：偏差法粒子 E 可正可负）
%   E_out : core  -> ghost 的总能量
%   E_net = E_in - E_out
%
% 依赖：attach_isothermal_ghosts_ 已把 mesh_g.core = [ix0,ix1; iy0,iy1; iz0,iz1] 写好，
%       且 ghost 只有“一层”：x 向 ghost 紧贴的层就是 ix0-1 或 ix1+1（y/z 同理）。
%
% 用法：
%   p_before = state.p;           % 飞行前拷贝
%   [state,~] = particle_fly_(...);   % 纯推进
%   iface = tally_iface_by_cid_(p_before, state.p, mesh_g, opts);

  if nargin<4, opts = struct(); end

  mk = @() struct('E_in',0,'E_out',0,'E_net',0,'N_in',0,'N_out',0);
  iface = struct('x_min',mk(),'x_max',mk(),'y_min',mk(), ...
                 'y_max',mk(),'z_min',mk(),'z_max',mk());

  if isempty(p_before), return; end

  % 网格与“core 矩形框”
  Nx = mesh_g.Nx; Ny = mesh_g.Ny; Nz = mesh_g.Nz;
  core = mesh_g.core;   % struct: ix0,ix1,iy0,iy1,iz0,iz1
  ix0 = core.ix0; ix1 = core.ix1;
  iy0 = core.iy0; iy1 = core.iy1;
  iz0 = core.iz0; iz1 = core.iz1;

  % 每个方向上的“ghost 层”下标（若该面未加 ghost，则置为空）
  ix_gL = []; ix_gR = [];
  iy_gB = []; iy_gT = [];
  iz_gD = []; iz_gU = [];

  if ix0 > 1,   ix_gL = ix0 - 1; end   % x_min
  if ix1 < Nx,  ix_gR = ix1 + 1; end   % x_max
  if iy0 > 1,   iy_gB = iy0 - 1; end   % y_min
  if iy1 < Ny,  iy_gT = iy1 + 1; end   % y_max
  if iz0 > 1,   iz_gD = iz0 - 1; end   % z_min
  if iz1 < Nz,  iz_gU = iz1 + 1; end   % z_max

  % 映射：id -> before 粒子
  map_before = containers.Map('KeyType','uint32','ValueType','any');
  for i=1:numel(p_before)
      map_before(uint32(p_before(i).id)) = p_before(i);
  end

  % after 还活着的粒子：对比前后索引，判断是否穿过接口，并按面记账
  for i=1:numel(p_after)
      pa = p_after(i);
      pid = uint32(pa.id);
      if ~isKey(map_before, pid), continue; end  % 这一步不应有新生
      pb = map_before(pid);

      if ~(pb.cell>0), continue; end

      [ixb,iyb,izb] = ind2sub([Nx,Ny,Nz], pb.cell);
      [ixa,iya,iza] = ind2sub([Nx,Ny,Nz], pa.cell);

      % 是否在 core 框内
      in_core_b = (ixb>=ix0 && ixb<=ix1) && (iyb>=iy0 && iyb<=iy1) && (izb>=iz0 && izb<=iz1);
      in_core_a = (ixa>=ix0 && ixa<=ix1) && (iya>=iy0 && iya<=iy1) && (iza>=iz0 && iza<=iz1);
      if in_core_b == in_core_a
          % core→core 或 ghost→ghost，忽略
          continue;
      end

      % 粒子能量（带符号，取“出发端”的 E）
      E = get_particle_energy_(pb, opts);

      % —— ghost -> core：计 E_in ——（单面穿越假设）
      if ~in_core_b && in_core_a
          if ~isempty(ix_gL) && ixb==ix_gL && ixa==ix0
              iface.x_min = add_in_(iface.x_min, E);
          elseif ~isempty(ix_gR) && ixb==ix_gR && ixa==ix1
              iface.x_max = add_in_(iface.x_max, E);
          elseif ~isempty(iy_gB) && iyb==iy_gB && iya==iy0
              iface.y_min = add_in_(iface.y_min, E);
          elseif ~isempty(iy_gT) && iyb==iy_gT && iya==iy1
              iface.y_max = add_in_(iface.y_max, E);
          elseif ~isempty(iz_gD) && izb==iz_gD && iza==iz0
              iface.z_min = add_in_(iface.z_min, E);
          elseif ~isempty(iz_gU) && izb==iz_gU && iza==iz1
              iface.z_max = add_in_(iface.z_max, E);
          else
              % 不是紧贴的一层：说明 dt 太大或斜向跨越，保守忽略
          end
      % —— core -> ghost：计 E_out ——
      elseif in_core_b && ~in_core_a
          if ~isempty(ix_gL) && ixb==ix0 && ixa==ix_gL
              iface.x_min = add_out_(iface.x_min, E);
          elseif ~isempty(ix_gR) && ixb==ix1 && ixa==ix_gR
              iface.x_max = add_out_(iface.x_max, E);
          elseif ~isempty(iy_gB) && iyb==iy0 && iya==iy_gB
              iface.y_min = add_out_(iface.y_min, E);
          elseif ~isempty(iy_gT) && iyb==iy1 && iya==iy_gT
              iface.y_max = add_out_(iface.y_max, E);
          elseif ~isempty(iz_gD) && izb==iz0 && iza==iz_gD
              iface.z_min = add_out_(iface.z_min, E);
          elseif ~isempty(iz_gU) && izb==iz1 && iza==iz_gU
              iface.z_max = add_out_(iface.z_max, E);
          else
              % 同上，保守忽略
          end
      end
  end

  % 汇总净值
  tags = fieldnames(iface);
  for k=1:numel(tags)
      t = tags{k}; b = iface.(t);
      b.E_net = b.E_in - b.E_out;
      iface.(t) = b;
  end
end

% ---------- helpers ----------
function b = add_in_(b, E),  b.E_in  = b.E_in  + E; b.N_in  = b.N_in  + 1; end
function b = add_out_(b, E), b.E_out = b.E_out + E; b.N_out = b.N_out + 1; end

function E = get_particle_energy_(pp, opts)
  if isfield(pp,'E') && isfinite(pp.E), E = pp.E; return; end
  sgn = +1; if isfield(pp,'sgn') && isfinite(pp.sgn), sgn = sign(pp.sgn); if sgn==0, sgn=+1; end, end
  Eeff = 1e-18; if isfield(opts,'E_eff') && isfinite(opts.E_eff) && opts.E_eff>0, Eeff = opts.E_eff; end
  E = sgn * Eeff;
end

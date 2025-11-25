function [mesh_g, ghost] = attach_isothermal_ghosts_(mesh, opts)
% 根据 mesh.bc 的 isothermal 面，在外侧拼一层 ghost cell（厚度等于相邻器件第一格）。
% 扩展后网格记：
%   mesh_g.core = struct('ix0','ix1','iy0','iy1','iz0','iz1')  % 原器件区在扩展网格里的索引范围
%   mesh_g.bc   : 等温面的“新外侧面”统一改为 absorb；其他面沿用原 mesh.bc
%   ghost.ids   : 所有 ghost cell 的线性索引（扩展网格坐标系）
%   ghost.by_face.(x_min/x_max/...) : 分面 ghost cell 的线性索引
%
% 依赖：getfield_or（本文件底部提供）

  %#ok<*AGROW>
  if nargin<2, opts = struct(); end
  mesh_g = mesh; ghost = struct();

  X = mesh.x_edges(:); Y = mesh.y_edges(:); Z = mesh.z_edges(:);
  Nx = mesh.Nx; Ny = mesh.Ny; Nz = mesh.Nz;

  need_xmin = is_isothermal_(mesh.bc,'x_min');
  need_xmax = is_isothermal_(mesh.bc,'x_max');
  need_ymin = is_isothermal_(mesh.bc,'y_min');
  need_ymax = is_isothermal_(mesh.bc,'y_max');
  need_zmin = is_isothermal_(mesh.bc,'z_min');
  need_zmax = is_isothermal_(mesh.bc,'z_max');

  % —— X 方向 ——（先左后右，便于计算 core 偏移）
  if need_xmin
    dx1 = X(2)-X(1);
    mesh_g.x_edges = [X(1)-dx1; X];
    mesh_g.Nx = Nx+1;
  else
    mesh_g.x_edges = X;
    mesh_g.Nx = Nx;
  end
  if need_xmax
    dxN = X(end)-X(end-1);
    mesh_g.x_edges = [mesh_g.x_edges; X(end)+dxN];
    mesh_g.Nx = mesh_g.Nx + 1;
  end

  % —— Y 方向 ——
  if need_ymin
    dy1 = Y(2)-Y(1);
    mesh_g.y_edges = [Y(1)-dy1; Y];
    mesh_g.Ny = Ny+1;
  else
    mesh_g.y_edges = Y;
    mesh_g.Ny = Ny;
  end
  if need_ymax
    dyN = Y(end)-Y(end-1);
    mesh_g.y_edges = [mesh_g.y_edges; Y(end)+dyN];
    mesh_g.Ny = mesh_g.Ny + 1;
  end

  % —— Z 方向 ——
  if need_zmin
    dz1 = Z(2)-Z(1);
    mesh_g.z_edges = [Z(1)-dz1; Z];
    mesh_g.Nz = Nz+1;
  else
    mesh_g.z_edges = Z;
    mesh_g.Nz = Nz;
  end
  if need_zmax
    dzN = Z(end)-Z(end-1);
    mesh_g.z_edges = [mesh_g.z_edges; Z(end)+dzN];
    mesh_g.Nz = mesh_g.Nz + 1;
  end

  % ========= 写 core 索引（原 Nx×Ny×Nz 在扩展网格中的包络）=========
  % 注意：若在“最小侧”加了 ghost，那么核心起点索引 +1
  ix0 = 1 + double(need_xmin);
  iy0 = 1 + double(need_ymin);
  iz0 = 1 + double(need_zmin);
  ix1 = ix0 + Nx - 1;
  iy1 = iy0 + Ny - 1;
  iz1 = iz0 + Nz - 1;
  mesh_g.core = struct('ix0',ix0,'ix1',ix1,'iy0',iy0,'iy1',iy1,'iz0',iz0,'iz1',iz1);

  % ========= 更新边界条件 =========
  % 原器件的等温面现在变成“内部接口”；真正的新外侧面统一设为 absorb。
  bc_new = mesh.bc;               % 先拷贝
  if need_xmin, bc_new.x_min = struct('type','absorb'); end
  if need_xmax, bc_new.x_max = struct('type','absorb'); end
  if need_ymin, bc_new.y_min = struct('type','absorb'); end
  if need_ymax, bc_new.y_max = struct('type','absorb'); end
  if need_zmin, bc_new.z_min = struct('type','absorb'); end
  if need_zmax, bc_new.z_max = struct('type','absorb'); end
  mesh_g.bc = bc_new;

  % ========= 生成 ghost cell 索引表 =========
  ghost = index_ghost_cells_(Nx,Ny,Nz, mesh_g, ...
            [need_xmin need_xmax need_ymin need_ymax need_zmin need_zmax], ...
            mesh_g.core);
end

% ---------------------------- 工具们 ----------------------------
function tf = is_isothermal_(bc, tag)
  tf = isfield(bc,tag) && isstruct(bc.(tag)) && isfield(bc.(tag),'type') ...
       && strcmpi(bc.(tag).type,'isothermal');
end

function ghost = index_ghost_cells_(Nx,Ny,Nz, mesh_g, needs, core)
% 在“扩展网格坐标系”下返回各面的 ghost cell 线性下标
  NGx = mesh_g.Nx; NGy = mesh_g.Ny; NGz = mesh_g.Nz;

  need_xmin=needs(1); need_xmax=needs(2);
  need_ymin=needs(3); need_ymax=needs(4);
  need_zmin=needs(5); need_zmax=needs(6);

  ghost = struct('ids',[],'by_face',struct(), 'flags', struct( ...
      'x_min',need_xmin,'x_max',need_xmax,'y_min',need_ymin,'y_max',need_ymax,'z_min',need_zmin,'z_max',need_zmax));

  add_ids = @(I) sub2ind([NGx,NGy,NGz], I{:});

  % 便于阅读
  ix0=core.ix0; ix1=core.ix1; iy0=core.iy0; iy1=core.iy1; iz0=core.iz0; iz1=core.iz1;

  if need_xmin
    % x=ix0-1 是 ghost 层；覆盖核心在 y,z 的投影
    I = { (ix0-1) * ones((iy1-iy0+1)*(iz1-iz0+1),1), ...
          repelem( (iy0:iy1).', (iz1-iz0+1) ), ...
          repmat( (iz0:iz1).', (iy1-iy0+1), 1 ) };
    ghost.by_face.x_min = add_ids(I);
  end
  if need_xmax
    % x=ix1+1
    I = { (ix1+1) * ones((iy1-iy0+1)*(iz1-iz0+1),1), ...
          repelem( (iy0:iy1).', (iz1-iz0+1) ), ...
          repmat( (iz0:iz1).', (iy1-iy0+1), 1 ) };
    ghost.by_face.x_max = add_ids(I);
  end

  if need_ymin
    % y=iy0-1
    I = { repelem( (ix0:ix1).', (iz1-iz0+1) ), ...
          (iy0-1) * ones((ix1-ix0+1)*(iz1-iz0+1),1), ...
          repmat( (iz0:iz1).', (ix1-ix0+1), 1 ) };
    ghost.by_face.y_min = add_ids(I);
  end
  if need_ymax
    % y=iy1+1
    I = { repelem( (ix0:ix1).', (iz1-iz0+1) ), ...
          (iy1+1) * ones((ix1-ix0+1)*(iz1-iz0+1),1), ...
          repmat( (iz0:iz1).', (ix1-ix0+1), 1 ) };
    ghost.by_face.y_max = add_ids(I);
  end

  if need_zmin
    % z=iz0-1
    I = { repelem( (ix0:ix1).', (iy1-iy0+1) ), ...
          repmat( (iy0:iy1).', (ix1-ix0+1), 1 ), ...
          (iz0-1) * ones((ix1-ix0+1)*(iy1-iy0+1),1) };
    ghost.by_face.z_min = add_ids(I);
  end
  if need_zmax
    % z=iz1+1
    I = { repelem( (ix0:ix1).', (iy1-iy0+1) ), ...
          repmat( (iy0:iy1).', (ix1-ix0+1), 1 ), ...
          (iz1+1) * ones((ix1-ix0+1)*(iy1-iy0+1),1) };
    ghost.by_face.z_max = add_ids(I);
  end

  % 汇总所有 ghost cell
  ghost.ids = [];
  fns = fieldnames(ghost.by_face);
  for k=1:numel(fns)
    ghost.ids = [ghost.ids; ghost.by_face.(fns{k})(:)];
  end
  ghost.ids = unique(ghost.ids(:));
end

function v = getfield_or(s, name, default_v)
  if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
      v = s.(name);
  else
      v = default_v;
  end
end

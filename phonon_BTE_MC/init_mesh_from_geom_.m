function mesh = init_mesh_from_geom_(cs)
% init_mesh_from_geom_ 生成规则笛卡尔网格并提供 O(1) 单元定位器（快路径友好）
%
% INPUT (geom):
%   cs.geom.L   = [Lx Ly Lz]   尺寸(米)
%   cs.mesh.Nx, .Ny, .Nz       各向 cell 个数 (>=1)
%   cs.bc       边界条件结构（可选），原样挂到 mesh.bc
%
% OUTPUT (mesh):
%   基本：.L, .Nx,.Ny,.Nz,.Nc, .dx,.dy,.dz,.hmin, .centers(Nc×3), .vol(Nc×1), .boxes(Nc×6)
%   索引器：.to_id, .from_id, .point2id
%   额外（用于加速）：.domain_box(1×6), .x_edges(1×Nx+1), .y_edges, .z_edges
%   （保持向后兼容，其它模块无需修改）

  % --------- 校验与基础量 ---------
  assert(isfield(cs.geom,'L') && numel(cs.geom.L)==3, 'geom.L = [Lx Ly Lz] 必须给出');
  Lx = cs.geom.L(1); Ly = cs.geom.L(2); Lz = cs.geom.L(3);
  assert(all([Lx,Ly,Lz] > 0), '几何尺寸必须为正');

  reqFields = {'Nx','Ny','Nz'};
  for f = reqFields
    assert(isfield(cs.mesh,f{1}), '缺少 geom.%s', f{1});
    assert(cs.mesh.(f{1})>=1 && mod(cs.mesh.(f{1}),1)==0, 'geom.%s 必须为正整数', f{1});
  end
  Nx = cs.mesh.Nx; Ny = cs.mesh.Ny; Nz = cs.mesh.Nz;
  Nc = Nx*Ny*Nz;

  if isfield(cs.mesh,'x_nodes') && isfield(cs.mesh,'y_nodes') && isfield(cs.mesh,'z_nodes')
      x_edges = cs.mesh.x_nodes(:).';
      y_edges = cs.mesh.y_nodes(:).';
      z_edges = cs.mesh.z_nodes(:).';
  else
      x_edges = linspace(0, Lx, Nx+1);
      y_edges = linspace(0, Ly, Ny+1);
      z_edges = linspace(0, Lz, Nz+1);
  end

  dx_vec = diff(x_edges); dy_vec = diff(y_edges); dz_vec = diff(z_edges);
  hmin = min([dx_vec(:); dy_vec(:); dz_vec(:)]);

  xc = 0.5*(x_edges(1:end-1) + x_edges(2:end));
  yc = 0.5*(y_edges(1:end-1) + y_edges(2:end));
  zc = 0.5*(z_edges(1:end-1) + z_edges(2:end));
  [Yc,Xc,Zc] = ndgrid(yc, xc, zc);
  centers = [Xc(:) Yc(:) Zc(:)];

  [Xmin,Ymin,Zmin] = ndgrid(x_edges(1:end-1), y_edges(1:end-1), z_edges(1:end-1));
  [Xmax,Ymax,Zmax] = ndgrid(x_edges(2:end),   y_edges(2:end),   z_edges(2:end));
  boxes = [Xmin(:) Xmax(:) Ymin(:) Ymax(:) Zmin(:) Zmax(:)];

  [DX,DY,DZ] = ndgrid(dx_vec, dy_vec, dz_vec);
  vol = (DX(:).*DY(:).*DZ(:));

  % --------- 索引器（O(1)）---------
  to_id   = @(i,j,k) sub2ind([Nx,Ny,Nz], i, j, k);
  from_id = @(id) local_from_id_(id, [Nx,Ny,Nz]);
  if all(diff(diff(x_edges))==0) && all(diff(diff(y_edges))==0) && all(diff(diff(z_edges))==0)
      dx = dx_vec(1); dy = dy_vec(1); dz = dz_vec(1);
      point2id = @(x,y,z) local_point2id_regular_(x,y,z, Lx,Ly,Lz, dx,dy,dz, Nx,Ny,Nz, to_id);
  else
      point2id = @(x,y,z) local_point2id_rectilinear_(x,y,z, x_edges,y_edges,z_edges, Nx,Ny,Nz, to_id);
  end

  % --------- 打包 mesh ---------
  mesh = struct();
  mesh.L   = [Lx,Ly,Lz];
  mesh.Nx  = Nx; mesh.Ny = Ny; mesh.Nz = Nz; mesh.Nc = Nc;
  mesh.dx  = mean(dx_vec); mesh.dy = mean(dy_vec); mesh.dz = mean(dz_vec); mesh.hmin = hmin;

  mesh.centers = centers;      % (Nc×3)
  mesh.vol     = vol;          % (Nc×1)
  mesh.boxes   = boxes;        % (Nc×6)
  mesh.cell_vol= vol;

  mesh.to_id     = to_id;      % 函数句柄
  mesh.from_id   = from_id;    % 函数句柄
  mesh.point2id  = point2id;   % 函数句柄

  % ---- 关键：给快路径提供的额外字段 ----
  mesh.domain_box = [x_edges(1), x_edges(end), y_edges(1), y_edges(end), z_edges(1), z_edges(end)];
  mesh.x_edges    = x_edges;
  mesh.y_edges    = y_edges;
  mesh.z_edges    = z_edges;

  % 边界条件原样随 mesh 带走
  if isfield(cs,'bc'), mesh.bc = cs.bc; else, mesh.bc = struct(); end
end

% =================== 本地小工具 ===================

function ijk = local_from_id_(id, sz)
  [i,j,k] = ind2sub(sz, id);
  ijk = [i,j,k];
end

function id = local_point2id_regular_(x,y,z, Lx,Ly,Lz, dx,dy,dz, Nx,Ny,Nz, to_id)
  % clamp 到域内，然后 floor 定位到 (i,j,k)
  x = max(0, min(Lx, x));
  y = max(0, min(Ly, y));
  z = max(0, min(Lz, z));
  if x >= Lx, ii = Nx; else, ii = floor(x/dx) + 1; end
  if y >= Ly, jj = Ny; else, jj = floor(y/dy) + 1; end
  if z >= Lz, kk = Nz; else, kk = floor(z/dz) + 1; end
  ii = min(max(ii,1), Nx);
  jj = min(max(jj,1), Ny);
  kk = min(max(kk,1), Nz);
  id = to_id(ii,jj,kk);
end

function id = local_point2id_rectilinear_(x,y,z, X,Y,Z, Nx,Ny,Nz, to_id)
  if x < X(1) || x > X(end) || y < Y(1) || y > Y(end) || z < Z(1) || z > Z(end)
      x = min(max(x, X(1)), X(end));
      y = min(max(y, Y(1)), Y(end));
      z = min(max(z, Z(1)), Z(end));
  end
  ix = find(X <= x, 1, 'last');
  iy = find(Y <= y, 1, 'last');
  iz = find(Z <= z, 1, 'last');
  ix = min(max(ix,1), Nx);
  iy = min(max(iy,1), Ny);
  iz = min(max(iz,1), Nz);
  id = to_id(ix, iy, iz);
end

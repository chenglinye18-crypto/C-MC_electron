function state = reindex_particles_to_mesh_(state, mesh_g)
% 用扩展后的网格 mesh_g 重新计算所有粒子的 cell 索引
% 位置 (x,y,z) 不变，但会被夹紧到对应 cell 内部的 [edge+eps, edge_next-eps]

  if isempty(state.p), return; end

  X = mesh_g.x_edges(:); Y = mesh_g.y_edges(:); Z = mesh_g.z_edges(:);
  Nx = mesh_g.Nx; Ny = mesh_g.Ny; Nz = mesh_g.Nz;

  % 批量取坐标
  x = [state.p.x].';  y = [state.p.y].';  z = [state.p.z].';

  % 先把坐标夹在域内，避免 discretize 落到边外
  epsl = 1e-12;
  x = min(max(x, X(1)+epsl), X(end)-epsl);
  y = min(max(y, Y(1)+epsl), Y(end)-epsl);
  z = min(max(z, Z(1)+epsl), Z(end)-epsl);

  % 用 mesh_g 的 edges 重新离散化
  ix = max(1, min(discretize(x, X), Nx));
  iy = max(1, min(discretize(y, Y), Ny));
  iz = max(1, min(discretize(z, Z), Nz));
  cid_new = int32(sub2ind([Nx,Ny,Nz], ix, iy, iz));

  % 再把坐标轻微夹紧到对应 cell 内，避免刚好踩到面
  x = min(max(x, X(ix)+epsl), X(ix+1)-epsl);
  y = min(max(y, Y(iy)+epsl), Y(iy+1)-epsl);
  z = min(max(z, Z(iz)+epsl), Z(iz+1)-epsl);

  % 回写
  for k = 1:numel(state.p)
      state.p(k).x    = x(k);
      state.p(k).y    = y(k);
      state.p(k).z    = z(k);
      state.p(k).cell = cid_new(k);
  end
end

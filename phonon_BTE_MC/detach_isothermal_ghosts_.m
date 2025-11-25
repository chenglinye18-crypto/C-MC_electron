function [state, mesh] = detach_isothermal_ghosts_(state, mesh_g, ghost, mesh_orig)
% DETACH_ISOTHERMAL_GHOSTS_
% 1) 删除所有处于 ghost cells 的粒子
% 2) 将剩余粒子的 cell 重新映射为 mesh_orig 的编号体系（位置不变，夹紧到 cell 内）
% 3) 返回 mesh=mesh_orig

  % --- 空场景：直接恢复网格 ---
  if isempty(state.p)
    mesh = mesh_orig; 
    return; 
  end

  % --- 标记 ghost cells ---
  NG = mesh_g.Nx * mesh_g.Ny * mesh_g.Nz;
  is_ghost_cell = false(NG,1);
  if isfield(ghost,'ids') && ~isempty(ghost.ids)
      gids = ghost.ids(:);
      gids = gids(gids>=1 & gids<=NG);
      is_ghost_cell(gids) = true;
  end

  % --- 过滤：仅保留 core 粒子 ---
  keep = true(numel(state.p),1);
  for i = 1:numel(state.p)
      ci = state.p(i).cell;
      if ci>=1 && ci<=NG
          keep(i) = ~is_ghost_cell(ci);
      else
          % 落在非法 cell 的粒子，保守起见也保留，后面用位置重算 cid
          keep(i) = true;
      end
  end
  state.p = state.p(keep);

  % --- 用原网格把剩余粒子 re-index 到 mesh_orig 的 cell，并夹紧位置 ---
  state = reindex_particles_to_mesh_local_(state, mesh_orig);

  % --- 恢复网格 ---
  mesh = mesh_orig;
end

% ================== 本地小工具：用目标网格重算 cid ==================
function state = reindex_particles_to_mesh_local_(state, mesh_target)
  if isempty(state.p), return; end

  X = mesh_target.x_edges(:);
  Y = mesh_target.y_edges(:);
  Z = mesh_target.z_edges(:);
  Nx = mesh_target.Nx; 
  Ny = mesh_target.Ny; 
  Nz = mesh_target.Nz;

  % 批量取坐标
  x = [state.p.x].';  
  y = [state.p.y].';  
  z = [state.p.z].';

  % 夹到目标域内，避免落到边外
  epsl = 1e-12;
  x = min(max(x, X(1)+epsl), X(end)-epsl);
  y = min(max(y, Y(1)+epsl), Y(end)-epsl);
  z = min(max(z, Z(1)+epsl), Z(end)-epsl);

  % 离散化到 cell
  ix = max(1, min(discretize(x, X), Nx));
  iy = max(1, min(discretize(y, Y), Ny));
  iz = max(1, min(discretize(z, Z), Nz));
  cid_new = int32(sub2ind([Nx,Ny,Nz], ix, iy, iz));

  % 再把坐标轻微夹紧到 cell 内部，避免“踩面”问题
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

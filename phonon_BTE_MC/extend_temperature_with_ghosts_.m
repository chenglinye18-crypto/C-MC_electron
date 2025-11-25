function Tcell_g = extend_temperature_with_ghosts_(Tcell_core, mesh, mesh_g, ghost, bc, LUT)
% 将核心网格温度扩展到含 ghost 的网格
% 规则：
%   1) 恒温面 ghost 切片 = Tb（若给 U 且提供 LUT，则 Tb = LUT.inv(U)）
%   2) 非恒温面 ghost 切片 = 邻近内层一片（拷贝）
%   3) 角/棱相交：若多恒温面相交，取算术平均；若混合恒温+拷贝，恒温优先（覆盖拷贝）
%
% 返回：与 mesh_g 线性下标一致的列向量

  % 尺寸
  Nx  = mesh.Nx;   Ny  = mesh.Ny;   Nz  = mesh.Nz;
  Nxg = mesh_g.Nx; Nyg = mesh_g.Ny; Nzg = mesh_g.Nz;

  % pad 标志
  pad = getfield_or(ghost, 'pad', struct());
  gx_min = getfield_or(pad,'x_min',0);  gx_max = getfield_or(pad,'x_max',0);
  gy_min = getfield_or(pad,'y_min',0);  gy_max = getfield_or(pad,'y_max',0);
  gz_min = getfield_or(pad,'z_min',0);  gz_max = getfield_or(pad,'z_max',0);

  % reshape 核心温度并嵌入扩展 3D 数组
  Tc = reshape(Tcell_core, [Nx, Ny, Nz]);
  Tg_sum = nan(Nxg, Nyg, Nzg);      % 叠加值
  Tg_cnt = zeros(Nxg, Nyg, Nzg);    % 权重/计数

  ix0 = 1 + gx_min;   ix1 = gx_min + Nx;
  iy0 = 1 + gy_min;   iy1 = gy_min + Ny;
  iz0 = 1 + gz_min;   iz1 = gz_min + Nz;

  % core 放入
  Tg_sum(ix0:ix1, iy0:iy1, iz0:iz1) = Tc;
  Tg_cnt(ix0:ix1, iy0:iy1, iz0:iz1) = 1;

  % ---- 小工具：向切片“叠加”或者“覆盖”写入（覆盖用于恒温优先） ----
  function add_slice_x(i, val, overwrite)
      if overwrite
          Tg_sum(i,:,:) = 0; Tg_cnt(i,:,:) = 0;
      end
      % 标量自动扩展即可
      prev = Tg_sum(i,:,:);
      if any(isnan(prev(:))), prev(isnan(prev)) = 0; end
      Tg_sum(i,:,:) = prev + val;
      Tg_cnt(i,:,:) = Tg_cnt(i,:,:) + 1;
  end
  function add_slice_y(j, val, overwrite)
      if overwrite
          Tg_sum(:,j,:) = 0; Tg_cnt(:,j,:) = 0;
      end
      prev = Tg_sum(:,j,:);
      if any(isnan(prev(:))), prev(isnan(prev)) = 0; end
      Tg_sum(:,j,:) = prev + val;
      Tg_cnt(:,j,:) = Tg_cnt(:,j,:) + 1;
  end
  function add_slice_z(k, val, overwrite)
      if overwrite
          Tg_sum(:,:,k) = 0; Tg_cnt(:,:,k) = 0;
      end
      prev = Tg_sum(:,:,k);
      if any(isnan(prev(:))), prev(isnan(prev)) = 0; end
      Tg_sum(:,:,k) = prev + val;
      Tg_cnt(:,:,k) = Tg_cnt(:,:,k) + 1;
  end

  % ---- 取恒温温度（支持 Tb 或 U→T 反解） ----
  function Tb = face_Tb(face)
      if ~isfield(bc, face), Tb = nan; return; end
      F = bc.(face);
      if isfield(F,'T') && isfinite(F.T)
          Tb = F.T;
      elseif isfield(F,'U') && isfinite(F.U)
          if nargin>=6 && ~isempty(LUT) && isfield(LUT,'inv')
              Tb = LUT.inv(F.U);
          else
              error('extend_temperature_with_ghosts_: 面 %s 提供了 U 但未提供 LUT 反解 T。', face);
          end
      else
          Tb = nan;
      end
  end

  % ================= X 面 =================
  if gx_min == 1
      if is_iso_(bc,'x_min')
          add_slice_x(1, face_Tb('x_min'), true);     % 恒温覆盖
      else
          add_slice_x(1, squeeze(Tg_sum(2,:,:)), false);  % 拷贝内层
      end
  end
  if gx_max == 1
      if is_iso_(bc,'x_max')
          add_slice_x(Nxg, face_Tb('x_max'), true);
      else
          add_slice_x(Nxg, squeeze(Tg_sum(Nxg-1,:,:)), false);
      end
  end

  % ================= Y 面 =================
  if gy_min == 1
      if is_iso_(bc,'y_min')
          add_slice_y(1, face_Tb('y_min'), true);
      else
          add_slice_y(1, squeeze(Tg_sum(:,2,:)), false);
      end
  end
  if gy_max == 1
      if is_iso_(bc,'y_max')
          add_slice_y(Nyg, face_Tb('y_max'), true);
      else
          add_slice_y(Nyg, squeeze(Tg_sum(:,Nyg-1,:)), false);
      end
  end

  % ================= Z 面 =================
  if gz_min == 1
      if is_iso_(bc,'z_min')
          add_slice_z(1, face_Tb('z_min'), true);
      else
          add_slice_z(1, squeeze(Tg_sum(:,:,2)), false);
      end
  end
  if gz_max == 1
      if is_iso_(bc,'z_max')
          add_slice_z(Nzg, face_Tb('z_max'), true);
      else
          add_slice_z(Nzg, squeeze(Tg_sum(:,:,Nzg-1)), false);
      end
  end

  % 归一化（角/棱相交求平均；core 与拷贝处计数为 1）
  mask = Tg_cnt > 0;
  Tg = Tg_sum;
  Tg(mask) = Tg_sum(mask) ./ Tg_cnt(mask);

  % 极端兜底：若仍有 NaN（理论上不会），填 core 均值
  if any(isnan(Tg), 'all')
      Tg(isnan(Tg)) = mean(Tc(:));
  end

  Tcell_g = Tg(:);
end

% --------- 小工具 ---------
function v = getfield_or(s, name, default_v)
  if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
      v = s.(name);
  else
      v = default_v;
  end
end

function tf = is_iso_(bc, tag)
  tf = isfield(bc, tag) && isfield(bc.(tag),'type') && strcmpi(bc.(tag).type,'isothermal') ...
       && ( (isfield(bc.(tag),'T') && isfinite(bc.(tag).T)) || isfield(bc.(tag),'U') );
end

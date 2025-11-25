function heat = load_heat_from_finfet(case_dir, step, varargin)
% LOAD_HEAT_FROM_FINFET  Convert electron MC heat#### to volumetric sources.
%
%   heat = load_heat_from_finfet('finfet', 6999);
%
% Parameters:
%   case_dir : path to finfet directory (containing lgrid.txt & data/)
%   step     : integer (e.g., 6999) selecting data/heat<step>
%
% Optional name-value pairs:
%   'Dt'             : electron MC time step (default 1e-16 s)
%   'StatHeatStep'   : stat_heat_step value (default 5000)
%   'Include'        : 'both' | 'electron' | 'hole'
%   'Scale'          : coordinate scale (default 1e-6, µm->m)
%   'SignConvention' : +1 to keep original sign, -1 to flip (default -1)
%
% Returns struct:
%   heat.qvol       : Nc×1 vector of volumetric power density [W/m^3]
%   heat.qvol_grid  : (Nx×Ny×Nz) array
%   heat.edges      : {x_edges,y_edges,z_edges}
%   heat.total_power: integral power [W]
%   heat.meta       : additional info (dt, stat_heat_step, file,...)

  arguments
      case_dir (1,:) char
      step (1,1) double {mustBeInteger, mustBeNonnegative}
      varargin.Name
      varargin.Value
  end

  ip = inputParser;
  addParameter(ip, 'Dt', 1e-16);
  addParameter(ip, 'StatHeatStep', 5000);
  addParameter(ip, 'Include', 'both');
  addParameter(ip, 'Scale', 1e-6);
  addParameter(ip, 'SignConvention', -1);
  parse(ip, varargin{:});
  opt = ip.Results;

  heat_file = fullfile(case_dir, 'data', sprintf('heat%d', step));
  if ~isfile(heat_file)
      error('load_heat_from_finfet: 未找到 %s', heat_file);
  end

  lgrid_path = fullfile(case_dir, 'lgrid.txt');
  grid = read_lgrid_nodes(lgrid_path, 'Scale', opt.Scale);
  x_nodes = grid.x_m;
  y_nodes = grid.y_m;
  z_nodes = grid.z_m;

  Nx_nodes = numel(x_nodes);
  Ny_nodes = numel(y_nodes);
  Nz_nodes = numel(z_nodes);
  if any([Nx_nodes,Ny_nodes,Nz_nodes] < 2)
      error('load_heat_from_finfet: 节点数量不足。');
  end

  raw = readmatrix(heat_file);
  if size(raw,2) < 7
      error('load_heat_from_finfet: heat 文件列数不足 (需要 >=7)');
  end

  idx_i = raw(:,1) + 1; % convert to 1-based
  idx_j = raw(:,2) + 1;
  idx_k = raw(:,3) + 1;
  if max(idx_i) > Nx_nodes || max(idx_j) > Ny_nodes || max(idx_k) > Nz_nodes
      error('load_heat_from_finfet: heat 文件索引超出范围');
  end

  switch lower(opt.Include)
      case 'electron'
          q_raw = raw(:,5);
      case 'hole'
          q_raw = raw(:,7);
      otherwise
          q_raw = raw(:,5) + raw(:,7);
  end
  q_raw = opt.SignConvention * q_raw;

  % eV/cm^3 per stat window -> W/m^3
  J_per_eV = 1.602176634e-19;
  scale_density = 1e6;  % cm^-3 -> m^-3
  scale_time = opt.Dt * opt.StatHeatStep;
  q_scale = J_per_eV * scale_density / max(scale_time, eps);

  node_field = zeros(Nx_nodes, Ny_nodes, Nz_nodes);
  lin_indices = sub2ind([Nx_nodes, Ny_nodes, Nz_nodes], idx_i, idx_j, idx_k);
  node_field = accumarray(lin_indices, q_raw, [Nx_nodes*Ny_nodes*Nz_nodes, 1]);
  node_field = reshape(node_field, Nx_nodes, Ny_nodes, Nz_nodes);
  node_field = node_field * q_scale;

  % 节点 -> cell: 取 8 个节点平均
  cell_field = 0.125 * ( ...
      node_field(1:end-1, 1:end-1, 1:end-1) + ...
      node_field(2:end,   1:end-1, 1:end-1) + ...
      node_field(1:end-1, 2:end,   1:end-1) + ...
      node_field(1:end-1, 1:end-1, 2:end  ) + ...
      node_field(2:end,   2:end,   1:end-1) + ...
      node_field(2:end,   1:end-1, 2:end  ) + ...
      node_field(1:end-1, 2:end,   2:end  ) + ...
      node_field(2:end,   2:end,   2:end  ));

  dx = diff(x_nodes); dy = diff(y_nodes); dz = diff(z_nodes);
  [DX,DY,DZ] = ndgrid(dx, dy, dz);
  cell_vol = DX(:) .* DY(:) .* DZ(:);

  heat = struct();
  heat.qvol_grid = cell_field;
  heat.qvol = reshape(cell_field, [], 1);
  heat.cell_vol = cell_vol;
  heat.edges = {x_nodes, y_nodes, z_nodes};
  heat.total_power = sum(heat.qvol(:) .* cell_vol);
  heat.meta = struct('file', heat_file, ...
                     'dt', opt.Dt, ...
                     'stat_heat_step', opt.StatHeatStep, ...
                     'include', opt.Include);
end

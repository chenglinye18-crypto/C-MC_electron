function cs = setup_case_from_electron(lgrid_path, varargin)
% SETUP_CASE_FROM_ELECTRON  Build BTE geometry/mesh from finfet lgrid.txt.
%
%   cs = setup_case_from_electron('finfet/lgrid.txt', 'BC', bcStruct);
%
% Optional parameters:
%   'Scale'     : coordinate scale (default 1e-6, µm -> m)
%   'BC'        : struct with fields x_min/x_max/... (default all adiabatic)
%   'ShiftOrigin': logical, true -> shift coordinates so xmin=ymin=zmin=0

  p = inputParser;
  addParameter(p, 'Scale', 1e-6);
  addParameter(p, 'BC', []);
  addParameter(p, 'ShiftOrigin', true);
  parse(p, varargin{:});

  grid = read_lgrid_nodes(lgrid_path, 'Scale', p.Results.Scale);

  x_nodes = grid.x_m;
  y_nodes = grid.y_m;
  z_nodes = grid.z_m;
  if ~p.Results.ShiftOrigin
      x_nodes = grid.x_nodes * p.Results.Scale;
      y_nodes = grid.y_nodes * p.Results.Scale;
      z_nodes = grid.z_nodes * p.Results.Scale;
  end

  Nx = numel(x_nodes) - 1;
  Ny = numel(y_nodes) - 1;
  Nz = numel(z_nodes) - 1;
  if Nx <= 0 || Ny <= 0 || Nz <= 0
      error('setup_case_from_electron: 节点数量不足，无法构建网格');
  end

  geom.origin = [x_nodes(1), y_nodes(1), z_nodes(1)];
  geom.L = [x_nodes(end) - x_nodes(1), ...
            y_nodes(end) - y_nodes(1), ...
            z_nodes(end) - z_nodes(1)];
  geom.shape = 'box';

  mesh = struct();
  mesh.Nx = Nx; mesh.Ny = Ny; mesh.Nz = Nz;
  mesh.x_nodes = x_nodes;
  mesh.y_nodes = y_nodes;
  mesh.z_nodes = z_nodes;

  bc = p.Results.BC;
  if isempty(bc)
      faces = {'x_min','x_max','y_min','y_max','z_min','z_max'};
      for i=1:numel(faces)
          bc.(faces{i}).type = 'adiabatic';
      end
  end

  cs = struct();
  cs.units.length = 'm';
  cs.units.temp   = 'K';
  cs.geom = geom;
  cs.mesh = mesh;
  cs.bc   = bc;
end

function grid = read_lgrid_nodes(lgrid_path, varargin)
% READ_LGRID_NODES  Parse finfet/lgrid.txt to obtain node coordinates.
%
%   grid = read_lgrid_nodes('finfet/lgrid.txt');
%   grid.x_nodes, grid.y_nodes, grid.z_nodes  -- raw coordinates (µm)
%   grid.units = 'micron';
%
%   可选参数：
%       'Scale' : 线性尺度转换因子（默认 1e-6，将 µm 转为 m）
%
% 文件格式：
%   每个方向依次给出 “节点个数 + 节点坐标列表”，支持空行。

  p = inputParser;
  addParameter(p, 'Scale', 1e-6);
  parse(p, varargin{:});
  scale = p.Results.Scale;

  fid = fopen(lgrid_path, 'r');
  if fid < 0
      error('read_lgrid_nodes: 无法打开 %s', lgrid_path);
  end
  cleaner = onCleanup(@() fclose(fid));

  coords = cell(1,3);
  for dim = 1:3
      nNodes = local_next_integer(fid);
      values = zeros(nNodes,1);
      idx = 1;
      while idx <= nNodes
          line = fgetl(fid);
          if ~ischar(line)
              error('read_lgrid_nodes: 读取第 %d 个方向的节点时遇到 EOF', dim);
          end
          line = strtrim(line);
          if isempty(line)
              continue;
          end
          nums = sscanf(line, '%f');
          if isempty(nums)
              continue;
          end
          take = min(numel(nums), nNodes - idx + 1);
          values(idx:idx+take-1) = nums(1:take);
          idx = idx + take;
      end
      coords{dim} = values;
  end

  grid = struct();
  grid.units = 'micron';
  grid.x_nodes = coords{1};
  grid.y_nodes = coords{2};
  grid.z_nodes = coords{3};
  grid.scale = scale;
  grid.origin = [grid.x_nodes(1), grid.y_nodes(1), grid.z_nodes(1)] * scale;
  grid.x_m = (grid.x_nodes - grid.x_nodes(1)) * scale;
  grid.y_m = (grid.y_nodes - grid.y_nodes(1)) * scale;
  grid.z_m = (grid.z_nodes - grid.z_nodes(1)) * scale;
end

function val = local_next_integer(fid)
  while true
      line = fgetl(fid);
      if ~ischar(line)
          error('read_lgrid_nodes: 意外到达文件结尾');
      end
      line = strtrim(line);
      if isempty(line)
          continue;
      end
      nums = sscanf(line, '%f');
      if isempty(nums)
          continue;
      end
      val = round(nums(1));
      return;
  end
end

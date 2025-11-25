function cs = setup_case_cube_100nm(varargin)
% SETUP_CASE_CUBE_100NM  初始化一个 100 nm 立方体热输运案例（仅几何+边界条件）。
%
% 用法：
%   cs = setup_case_cube_100nm();            % 默认网格 32x32x32
%   cs = setup_case_cube_100nm('Nx',50);     % 自定义网格
%   cs = setup_case_cube_100nm('Nx',40,'Ny',30,'Nz',20);
%
% 输出结构体 cs，包含：
%   cs.units  — 字符串说明单位
%   cs.geom   — 几何定义（形状、尺寸）
%   cs.mesh   — 规则网格（dx,dy,dz、中心坐标等），便于后续离散
%   cs.bc     — 六个面的边界条件（等温/绝热/周期），温度单位 K
%
% 说明：
%   这里选 x=0 面为 400 K，x=Lx 面为 300 K，其余四面绝热。
%   若想对调冷热面，只需交换 bc.x_min.T 与 bc.x_max.T。

% ---------------- 参数解析（网格密度可改） ----------------
p = inputParser;
addParameter(p, 'Nx', 32);
addParameter(p, 'Ny', 32);
addParameter(p, 'Nz', 32);
parse(p, varargin{:});
Nx = p.Results.Nx; Ny = p.Results.Ny; Nz = p.Results.Nz;

% ---------------- 单位与几何（易于替换） ----------------
cs.units.length = 'm';
cs.units.temp   = 'K';

geom.shape   = 'box';                 % 形状类型：'box'（长方体/立方体）
geom.origin  = [0,0,0];               % 起点（左下前角）
geom.L       = [100,100,100]*1e-9;    % 尺寸 [m]：100 nm 立方体
% 以后换几何：
%   - 薄膜/矩形：仍用 'box'，只改 L
%   - 圆柱：geom.shape='cylinder'; 再添加半径/高度，扩展下方的 make_regular_mesh

cs.geom = geom;

% ---------------- 网格（规则网格，后续好离散） ----------------
cs.mesh = make_regular_mesh(geom, Nx, Ny, Nz);

% ---------------- 边界条件（清晰可读、易改） ----------------
% 约定：type = 'isothermal' | 'adiabatic' | 'periodic'
%       T    = 墙面温度（仅等温需要）
bc.x_min.type = 'isothermal';  bc.x_min.T = 350.0; % x=0   热端
bc.x_max.type = 'isothermal';  bc.x_max.T = 350.0; % x=Lx  冷端
%bc.y_min.type = 'isothermal';  bc.y_min.T = 401.0;
%bc.y_max.type = 'isothermal';  bc.y_max.T = 401.0;
%bc.z_min.type = 'isothermal';  bc.z_min.T = 401.0;
%bc.z_max.type = 'isothermal';  bc.z_max.T = 401.0;
%bc.x_min.type = 'adiabatic';
%bc.x_max.type = 'adiabatic';
bc.y_min.type = 'adiabatic';
bc.y_max.type = 'adiabatic';
bc.z_min.type = 'adiabatic';
bc.z_max.type = 'adiabatic';

cs.bc = bc;

% ---------------- 小结（便于检查） ----------------
print_case_summary(cs);

end

% ====== 子函数：规则网格（面心/体心都易拿） ======
function mesh = make_regular_mesh(geom, Nx, Ny, Nz)
assert(strcmpi(geom.shape,'box'), 'make_regular_mesh 目前只支持 box，其他形状可在此扩展。');
Lx = geom.L(1); Ly = geom.L(2); Lz = geom.L(3);

dx = Lx / Nx; dy = Ly / Ny; dz = Lz / Nz;

% 单元中心坐标（后续做有限体积/统计量更方便）
xc = linspace(dx/2, Lx - dx/2, Nx);
yc = linspace(dy/2, Ly - dy/2, Ny);
zc = linspace(dz/2, Lz - dz/2, Nz);
[Xc,Yc,Zc] = ndgrid(xc, yc, zc);

mesh.Nx = Nx; mesh.Ny = Ny; mesh.Nz = Nz;
mesh.dx = dx; mesh.dy = dy; mesh.dz = dz;
mesh.Lx = Lx; mesh.Ly = Ly; mesh.Lz = Lz;

mesh.xc = xc; mesh.yc = yc; mesh.zc = zc; % 1D 坐标
mesh.Xc = Xc; mesh.Yc = Yc; mesh.Zc = Zc; % 3D 网格（中心点）

% 六个面的面积（有时做通量/归一化会用到）
mesh.Ax = Ly * Lz;   % 与 x 正交的面
mesh.Ay = Lx * Lz;   % 与 y 正交的面
mesh.Az = Lx * Ly;   % 与 z 正交的面

% 体积
mesh.V  = Lx * Ly * Lz;
mesh.cellVolume = dx * dy * dz;
end

% ====== 子函数：打印小结，方便检查设置 ======
function print_case_summary(cs)
g = cs.geom; m = cs.mesh; bc = cs.bc;
fprintf('[geom] %s | L = (%.0f nm, %.0f nm, %.0f nm)\n', ...
  g.shape, g.L(1)*1e9, g.L(2)*1e9, g.L(3)*1e9);
fprintf('[mesh] Nx,Ny,Nz = %d, %d, %d | dx,dy,dz = (%.2g, %.2g, %.2g) m\n', ...
  m.Nx, m.Ny, m.Nz, m.dx, m.dy, m.dz);
fprintf('[bc] x_min: %s', bc.x_min.type);
if isfield(bc.x_min,'T'), fprintf(' (T=%.1f K)', bc.x_min.T); end
fprintf(' | x_max: %s', bc.x_max.type);
if isfield(bc.x_max,'T'), fprintf(' (T=%.1f K)\n', bc.x_max.T); else, fprintf('\n'); end
fprintf('[bc] y_min: %s', bc.y_min.type);
if isfield(bc.y_min,'T'), fprintf(' (T=%.1f K)', bc.y_min.T); end
fprintf(' | y_max: %s', bc.y_max.type);
if isfield(bc.y_max,'T'), fprintf(' (T=%.1f K)\n', bc.y_max.T); else, fprintf('\n'); end
fprintf('[bc] z_min: %s', bc.z_min.type);
if isfield(bc.z_min,'T'), fprintf(' (T=%.1f K)', bc.z_min.T); end
fprintf(' | z_max: %s', bc.z_max.type);
if isfield(bc.z_max,'T'), fprintf(' (T=%.1f K)\n', bc.z_max.T); else, fprintf('\n'); end

end

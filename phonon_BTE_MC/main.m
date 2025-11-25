clc;
clear all;
fprintf('*****************************************************\n');
fprintf('*          BTE Monte Carlo Simulator V1.0           *\n');
fprintf('*    Developed by Chenglin Ye, Peking University    *\n');
fprintf('*           Release Date: 20th Oct, 2025            *\n');
fprintf('*****************************************************\n');
cs = setup_case_cube_100nm('Nx',100,'Ny',1,'Nz',1);  % 得到几何+网格+边界
si   = mat_silicon_100();                               % 硅(100)色散
opts = mc_default_opts();                               
[Tp, p, out] = MC_solve_BTE(cs, si, opts);  
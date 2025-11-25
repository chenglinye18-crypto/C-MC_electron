function [Tprime, p, out] = MC_time_loop_BTE(mesh, spec, opts, state)
% MC_time_loop_BTE  能量包 MC-BTE 主循环（ghost-热库 + 飞行后接口能量核算）
% 依赖：
%   attach_isothermal_ghosts_, reindex_particles_to_mesh_, seed_ghost_particles_,
%   extend_temperature_with_ghosts_, precompute_relax_times_, step_pick_dt_,
%   particle_fly_, tally_iface_by_cid_, detach_isothermal_ghosts_,
%   particle_scattering_, update_temperature_from_energy_, build_E_T_lookup_

    % ---------------- 参数 ----------------
    Nc       = infer_Nc_(mesh);
    T0       = get_or(opts,'T0', 300);
    dt_min   = get_or(opts,'dt_min', 1e-15);
    dt_max   = get_or(opts,'dt_max', 1e-10);
    max_steps= get_or(opts,'max_steps', 5000);
    alpha_T  = max(min(get_or(opts,'T_underrelax', 1.0),1),0);
    scatter_on = get_or(opts,'scatter_on', true);

 %   %----------------- 热源 ---------------
 %   % ===== 体热源（保留模板，默认注释掉；以后需要再启用） =====
 %   % 网格边界
 %   % X = mesh.x_edges; Y = mesh.y_edges; Z = mesh.z_edges;
 %   % 取中间三分之一的盒子区域
 %   % ix1 = 1 + floor(mesh.Nx/3);    ix2 = mesh.Nx - floor(mesh.Nx/3);
 %   % iy1 = 1 + floor(mesh.Ny/3);    iy2 = mesh.Ny - floor(mesh.Ny/3);
 %   % iz1 = 1 + floor(mesh.Nz/3);    iz2 = mesh.Nz - floor(mesh.Nz/3);
 %   % box_bounds = [ X(ix1) X(ix2)  Y(iy1) Y(iy2)  Z(iz1) Z(iz2) ];  % [xmin xmax ymin ymax zmin zmax]
 %   % srcV = struct( ...
 %   %     'type',   'volume', ...
 %   %     'qvol',   1e14, ...                 % 体功率密度 [W/m^3]
 %   %     'E_eff',  opts.E_eff, ...           % 每个超粒子能量
 %   %     'region', struct('type','box','bounds',box_bounds) ...
 %   % );

    % ===== 面热源：位于 x 中间，覆盖整个 yz 面；向两侧同时散热 =====
    source_cfg = getfield_or(opts,'source', struct());
    use_volume_map = isfield(source_cfg,'qvol');
    if ~use_volume_map
        X = mesh.x_edges; Y = mesh.y_edges; Z = mesh.z_edges;
        x_mid   = 0.5*(X(1) + X(end));
        y_min   = Y(1);     y_max = Y(end);
        z_min   = Z(1);     z_max = Z(end);

        bounds_midYZ = [z_min z_max  y_min y_max  x_mid];

        qsurf_val = 1e10;   % [W/m^2]
        srcS_posX = struct( ...
            'type',  'surface', ...
            'qsurf', qsurf_val, ...
            'E_eff', [], ...
            'face',  struct('normal',[+1 0 0], 'bounds', bounds_midYZ) ...
        );
        srcS_negX = struct( ...
            'type',  'surface', ...
            'qsurf', qsurf_val, ...
            'E_eff', [], ...
            'face',  struct('normal',[-1 0 0], 'bounds', bounds_midYZ) ...
        );
    end


    % 日志
    logcfg       = get_or(opts,'log', struct());
    verbose_faces = get_or(logcfg,'verbose_faces', true);  % 是否打印逐面统计
    log_on       = get_or(logcfg,'on', true);
    print_every  = get_or(logcfg,'print_every', 1);
    header_every = get_or(logcfg,'header_every', 50);
    to_file      = get_or(logcfg,'to_file', false);
    logfile      = get_or(logcfg,'filename','mc_log.txt');
    if to_file, fid = fopen(logfile,'a'); if fid<0, fid=1; end, else, fid=1; end

    % ---- 收敛配置（可从 opts.conv 覆盖）----
    conv = struct();
    conv.min_steps = get_or(getfield_or(opts,'conv',struct()), 'min_steps', 10);     % 最少步数
    conv.n_consec  = get_or(getfield_or(opts,'conv',struct()), 'n_consec',  3);      % 连续通过轮数
    conv.tol_inf   = get_or(getfield_or(opts,'conv',struct()), 'tol_inf',   5e-2);   % K
    conv.tol_l2    = get_or(getfield_or(opts,'conv',struct()), 'tol_l2',    5e-2);   % K
    conv.tol_Enet  = get_or(getfield_or(opts,'conv',struct()), 'tol_Enet',  2e-18);  % J

    consec_ok = 0;   % 连续通过计数器


    hdr = [...
      '  step |      dt[s]   |  dT_inf[K] |  dT_L2[K] |    E_net[J] |' ...
      '   T_min[K] |  T_mean[K] |   T_max[K] |  pscat_max |      Np' newline];

    fmt = '%6d | %1.4e | %1.4e | %1.4e | %+.3e | %9.3f | %10.3f | %9.3f | %10.3f | %7d\n';

    % ---------------- 初始化 ----------------
    Tstar  = T0*ones(Nc,1);
    Tprime = Tstar;

    out = struct('dt_hist',[], 'T_inf_hist',[], 'T_l2_hist',[], 'pscat_max_hist',[], ...
                 'E_net_hist',[], 'dU_cells_hist',[], 'dU_alive_hist',[], 'resid_hist',[], ...
                 'iface_hist',[], 'nsteps',0, 'converged',false, 'Temperature_hist',[]);

    % 温度→能量密度查表（你的版本：LUT.T, LUT.U, LUT.inv）
    LUT = build_E_T_lookup_(spec, struct('T_min',1,'T_max',2000,'nT',2001));
    U_of_T = @(T) interp1(LUT.T, LUT.U, clamp_vec(T, LUT.T(1), LUT.T(end)), 'pchip');

    LUTq = build_q_T_lookup_(spec, struct('T_min',1,'T_max',2000,'nT',2001));

    % cell 体积（能量账本）
    Vc = cell_volumes_(mesh);

    % 初始能量
    U_alive_prev = particles_total_energy_(state, opts);       % 粒子（偏差：净偏差能量）
    U_cells_prev = sum( U_of_T(Tstar) .* Vc );                 % cell 绝对能量

    if log_on
        fprintf(fid,'[%s] MC BTE start. Ncells=%d, T0=%.2f K\n', datestr(now,'HH:MM:SS'), Nc, T0);
        fprintf(fid, hdr);
    end

    % ================== 时间循环 ==================
    for step = 1:max_steps
        draw_iter_snapshot(state, spec, mesh, step, opts);
        % ---- (0) ghost 热库（发射内化）----
        [mesh_g, ghost] = attach_isothermal_ghosts_(mesh, opts);
        state = reindex_particles_to_mesh_(state, mesh_g);
        state = seed_ghost_particles_(state, mesh_g, ghost, spec, opts, mesh.bc);
        %Tstar_g = extend_temperature_with_ghosts_(Tstar, mesh, mesh_g, ghost, mesh.bc, LUT);
        [Tstar_g, ~] = update_temperature_from_energy_(state, mesh_g, spec, opts, LUT);

        % ---- (A) 弛豫时间 ----
        r_tau = precompute_relax_times_(state, Tstar_g, opts, spec);
        if ~scatter_on, r_tau(:,:) = 0; end

        % ---- (B) 选 dt ----
        [dt, info_dt] = step_pick_dt_(mesh_g, spec, opts, state, r_tau);
        if ~isfinite(dt), dt = dt_min; end
        dt = min(max(dt, dt_min), dt_max);
        out.dt_hist(end+1,1) = dt;

        if use_volume_map
            newpV = spawn_volume_map_sources_(source_cfg, opts, mesh, spec, state, Tprime, LUT, LUTq, dt);
            if ~isempty(newpV)
                state.p = [state.p; newpV];
            end
        else
            srcS_posX.E_eff = state.WE;
            srcS_negX.E_eff = state.WE;
            newp1 = spawn_heat_source(opts, mesh, spec, state, Tprime, LUT, LUTq, srcS_posX, dt);
            state.p = [state.p; newp1];
            newp2 = spawn_heat_source(opts, mesh, spec, state, Tprime, LUT, LUTq, srcS_negX, dt);
            state.p = [state.p; newp2];
        end


        % ---- (C) 飞行 + 接口能量核算 ----
        p_before = state.p;
        [state, ~] = particle_fly_(state, mesh_g, dt, opts, spec);
        p_after = state.p;
        iface = tally_iface_by_cid_(p_before, p_after, mesh_g, opts);

        faces = {'x_min','x_max','y_min','y_max','z_min','z_max'};
        E_net_total = 0;
        for i=1:numel(faces)
            f = faces{i};
            if isfield(iface,f)
                E_net_total = E_net_total + iface.(f).E_net;
            end
        end
        out.E_net_hist(end+1,1) = E_net_total;
        out.iface_hist{end+1} = iface; 

        % ---- (D) 拆 ghost ----
        [state, mesh] = detach_isothermal_ghosts_(state, mesh_g, ghost, mesh);

        % ---- (E) 体散射 ----
        if scatter_on
            r_tau = precompute_relax_times_(state, Tstar, opts, spec);
            try
                [state, ~] = particle_scattering_(state, mesh, spec, opts, dt, Tstar, r_tau);
            catch
                state = particle_scattering_(state, mesh, spec, opts, dt, Tstar, r_tau);
            end
        end

        % ---- (F) 反推温度 ----
        [Tnew, ~] = update_temperature_from_energy_(state, mesh, spec, opts, LUT);
        Tprime = Tnew;
        out.Temperature_hist(end+1,:) = Tprime;

        % ---- (G) 能量账本 + 温差 ----
        U_alive_now = particles_total_energy_(state, opts);
        U_cells_now = sum( U_of_T(Tprime) .* Vc );

        dU_cells = U_cells_now - U_cells_prev;
        dU_alive = U_alive_now - U_alive_prev;
        resid    = dU_cells - dU_alive;

        out.dU_cells_hist(end+1,1) = dU_cells;
        out.dU_alive_hist(end+1,1) = dU_alive;
        out.resid_hist(end+1,1)    = resid;

        dT    = Tprime - Tstar;
        T_inf = norm(dT, inf);
        T_l2  = norm(dT)/sqrt(max(Nc,1));
        out.T_inf_hist(end+1,1) = T_inf;
        out.T_l2_hist(end+1,1)  = T_l2;

    
        pscat_max = NaN;
        if isstruct(info_dt) && isfield(info_dt,'p_scat_max') && ~isempty(info_dt.p_scat_max)
            pscat_max = info_dt.p_scat_max;
        end
        out.pscat_max_hist(end+1,1) = pscat_max;

        if log_on && (step==1 || mod(step, print_every)==0)
            % 每次输出数据行前先打一遍表头
            fprintf(fid, hdr);

            Tmin = min(Tprime); Tmean = mean(Tprime); Tmax = max(Tprime);
            fprintf(fid, fmt, step, dt, T_inf, T_l2, E_net_total, Tmin, Tmean, Tmax, ...
                    safe_num(pscat_max,0), numel(state.p));

            % 逐面统计（可关）
            if verbose_faces
                faces = {'x_min','x_max','y_min','y_max','z_min','z_max'};
                for ii = 1:numel(faces)
                    f = faces{ii};
                    if isfield(iface,f)
                        b = iface.(f);
                        fprintf(fid,'   [face %6s] E_in=%+.3e  E_out=%+.3e  E_net=%+.3e  (N_in=%d, N_out=%d)\n', ...
                                f, b.E_in, b.E_out, b.E_net, b.N_in, b.N_out);
                    end
                end
            end

            % 能量账本 & dt 诊断
            fprintf(fid,'   [ledger] dUc=%+.3e  dUa=%+.3e  resid=%+.3e  (J)\n', dU_cells, dU_alive, resid);
            if isstruct(info_dt) && (isfield(info_dt,'dt_cfl') || isfield(info_dt,'dt_prob'))
                fprintf(fid,'   [dt] dt_cfl=%1.3e  dt_prob=%1.3e\n', ...
                    safe_num(getfield_or(info_dt,'dt_cfl',NaN),NaN), ...
                    safe_num(getfield_or(info_dt,'dt_prob',NaN),NaN));
            end

            % 刷新输出（文件/控制台）
            if to_file && fid>1, fflush(fid); end
            drawnow limitrate;
        end


        % 欠松弛与能量基线更新
        if alpha_T < 1, Tstar = (1-alpha_T)*Tstar + alpha_T*Tprime; else, Tstar = Tprime; end
        U_cells_prev = U_cells_now;
        U_alive_prev = U_alive_now;

        out.nsteps = step;

          % ---- 收敛判定（最小步数 + 连续通过）----
        pass_now = (T_inf <= conv.tol_inf) && (T_l2 <= conv.tol_l2) ...
                   && (abs(E_net_total) <= conv.tol_Enet);

        if pass_now
            consec_ok = consec_ok + 1;
        else
            consec_ok = 0;
        end

        if (step >= conv.min_steps) && (consec_ok >= conv.n_consec)
            out.converged = true;
            out.nsteps    = step;
            if log_on
                fprintf(fid,'[%s] Converged at step %d: dT_inf=%1.3e, dT_L2=%1.3e, E_net=%+.3e\n', ...
                        datestr(now,'HH:MM:SS'), step, T_inf, T_l2, E_net_total);
            end
            break;
        end
    end 

    if log_on && to_file && fid>1, fclose(fid); end
    p = state.p;
end

% ======================= 工具函数（本文件内部） =======================
function p = spawn_volume_map_sources_(src_cfg, opts, mesh, spec, state, Tprime, LUT, LUTq, dt)
    qvol = src_cfg.qvol(:);
    if isfield(src_cfg,'mask') && ~isempty(src_cfg.mask)
        mask = logical(src_cfg.mask(:));
    else
        mask = abs(qvol) > 0;
    end
    cell_ids = find(mask);
    if isempty(cell_ids), p = struct([]); return; end
    E_eff = getfield_or(src_cfg,'E_eff', state.WE);
    p = struct([]);
    for idx = reshape(cell_ids,1,[])
        src_local = struct( ...
            'type',  'volume', ...
            'qvol',  qvol(idx), ...
            'E_eff', E_eff, ...
            'region', struct('type','cells','id', idx) ...
        );
        chunk = spawn_heat_source(opts, mesh, spec, state, Tprime, LUT, LUTq, src_local, dt);
        if ~isempty(chunk)
            if isempty(p)
                p = chunk;
            else
                p = [p; chunk]; %#ok<AGROW>
            end
        end
    end
end

function Nc = infer_Nc_(mesh)
    if isfield(mesh,'Nc') && ~isempty(mesh.Nc), Nc = mesh.Nc; return; end
    if all(isfield(mesh,{'Nx','Ny','Nz'})), Nc = mesh.Nx*mesh.Ny*mesh.Nz; return; end
    if isfield(mesh,'boxes'), Nc = size(mesh.boxes,1); return; end
    error('infer_Nc_: cannot infer Nc from mesh.');
end

function v = get_or(s, name, default_v)
    if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = default_v; end
end

function v = getfield_or(s, name, default_v)
    if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = default_v; end
end

function x = clamp_vec(x, a, b), x = min(max(x,a), b); end

function Vc = cell_volumes_(mesh)
    if all(isfield(mesh,{'Nx','Ny','Nz','x_edges','y_edges','z_edges'}))
        dx = diff(mesh.x_edges(:)); dy = diff(mesh.y_edges(:)); dz = diff(mesh.z_edges(:));
        Vc3 = reshape(dx,[],1).*reshape(dy,1,[]).*reshape(dz,1,1,[]);
        Vc  = Vc3(:);
    elseif isfield(mesh,'boxes') && ~isempty(mesh.boxes)
        bx = mesh.boxes; Vc = (bx(:,2)-bx(:,1)).*(bx(:,4)-bx(:,3)).*(bx(:,6)-bx(:,5));
    elseif isfield(mesh,'cell_vol') && ~isempty(mesh.cell_vol)
        Vc = mesh.cell_vol(:);
    else
        error('cell_volumes_: mesh lacks edges/boxes/cell_vol.');
    end
end

function E = particles_total_energy_(state, opts)
    if isempty(state.p), E = 0; return; end
    if isfield(state.p,'E') && ~isempty([state.p.E])
        E = sum([state.p.E]);
    else
        E_eff = get_or(opts,'E_eff',1e-18);
        E = E_eff * numel(state.p);
    end
end

function v = safe_num(v, fallback)
    if ~isfinite(v), v = fallback; end
end

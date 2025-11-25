function [is_ok, stat] = check_convergence_(hist, thresh)
% CHECK_CONVERGENCE_
% 兼容多种历史量；缺少的字段自动跳过相应判据。
%
% 输入 hist 可包含（任意子集）：
%   hist.T_inf(k)       % 步 k 的 ||ΔT||_inf
%   hist.T_l2(k)        % 步 k 的 ||ΔT||_2 / sqrt(Nc)
%   hist.E_emit(k)      % 步 k 发射能量 (J)
%   hist.E_absorb(k)    % 步 k 吸收能量 (J)
%   hist.U_alive(k)     % 步末域内总能量 (J)
%   hist.dt(k), hist.step(k)  % 可选
%
% 阈值 thresh（给哪个用哪个，未给则跳过该判据）：
%   .tol_inf   (K)        温度无穷范数阈值（窗口内最大值）
%   .tol_l2    (K)        温度 L2 阈值（窗口内最大值）
%   .tol_Ebal  (J/step)   能量收支残差均值阈值（窗口内）
%   .tol_dU    (J/step)   域能变化均值阈值（窗口内）
%   .window    (int,5)    统计窗口长度
%   .min_steps (int,5)    至少迭代步数才开始判定
%   .consecutive (int,1)  需要连续满足的次数（本函数仅返回“本次是否 OK”，连续计数放主循环）

  % ---- 参数 ----
  window      = getd(thresh,'window',5);
  min_steps   = getd(thresh,'min_steps',5);
  tol_inf     = getd(thresh,'tol_inf',[]);
  tol_l2      = getd(thresh,'tol_l2',[]);
  tol_Ebal    = getd(thresh,'tol_Ebal',[]);
  tol_dU      = getd(thresh,'tol_dU',[]);

  k = maxlen(hist);
  if k < min_steps
      is_ok = false;
      stat = base_stat(false);
      return;
  end
  i0  = max(1, k-window+1);
  idx = i0:k;

  ok_list = []; details = struct();

  % ---- 温度判据（用 T_inf / T_l2 历史，而非 T_curr/T_prev）----
  if isfield(hist,'T_inf') && ~isempty(tol_inf) && ~isempty(hist.T_inf)
      x = hist.T_inf(idx);
      ok_Tinf = all(isfinite(x)) && max(x) <= tol_inf;
      ok_list(end+1) = ok_Tinf; %#ok<AGROW>
      details.max_Tinf = max(x);
  end

  if isfield(hist,'T_l2') && ~isempty(tol_l2) && ~isempty(hist.T_l2)
      x = hist.T_l2(idx);
      ok_Tl2 = all(isfinite(x)) && max(x) <= tol_l2;
      ok_list(end+1) = ok_Tl2; %#ok<AGROW>
      details.max_Tl2 = max(x);
  end

  % ---- 边界能量平衡：E_emit - E_absorb ≈ ΔU ----
  if isfield(hist,'E_emit') && isfield(hist,'E_absorb') && ~isempty(tol_Ebal) ...
     && ~isempty(hist.E_emit) && ~isempty(hist.E_absorb)
      Eem = hist.E_emit(idx);
      Eab = hist.E_absorb(idx);
      if isfield(hist,'U_alive') && numel(hist.U_alive) >= k
          U  = hist.U_alive;
          % 与 dU 对齐：取窗口前一刻作为上一状态
          Useg = U(max(1,i0-1):k);
          if numel(Useg) >= 2
              dU = diff(Useg);
              dU = dU(end-numel(Eem)+1:end); % 与窗口长度对齐
          else
              dU = zeros(size(Eem));
          end
      else
          % 没有 U_alive 就只看 (E_emit - E_absorb) 的绝对量
          dU = zeros(size(Eem));
      end
      resid = Eem - Eab - dU;             % 能量残差
      mean_abs_resid = mean(abs(resid));
      ok_Ebal = isfinite(mean_abs_resid) && mean_abs_resid <= tol_Ebal;
      ok_list(end+1) = ok_Ebal; %#ok<AGROW>
      details.mean_abs_Eresid = mean_abs_resid;
  end

  % ---- 域能稳定：|ΔU| 的平均值小 ----
  if isfield(hist,'U_alive') && ~isempty(tol_dU) && numel(hist.U_alive) >= 2
      Useg = hist.U_alive(max(1,i0-1):k);
      dU   = diff(Useg);
      if isempty(dU), dU = 0; end
      mean_abs_dU = mean(abs(dU));
      ok_dU = isfinite(mean_abs_dU) && mean_abs_dU <= tol_dU;
      ok_list(end+1) = ok_dU; %#ok<AGROW>
      details.mean_abs_dU = mean_abs_dU;
  end

  % ---- 聚合 ----
  if isempty(ok_list)
      is_ok = false;
      stat = base_stat(false);
      return;
  end
  is_ok = all(ok_list);
  stat  = base_stat(is_ok);
  stat.details = details;
end

% ---------- helpers ----------
function v = getd(s, name, def)
  if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end

function L = maxlen(hist)
  L = 0;
  fns = fieldnames(hist);
  for i = 1:numel(fns)
      v = hist.(fns{i});
      if isnumeric(v), L = max(L, numel(v)); end
  end
end

function s = base_stat(ok)
  s = struct('ok_now', logical(ok), 'details', struct());
end

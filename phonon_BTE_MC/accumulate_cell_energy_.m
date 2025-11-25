function [E_cell, stats] = accumulate_cell_energy_(state, mesh)
% 汇总每个 cell 的总能量（J）
  Nc = size(mesh.boxes,1);
  cid = [state.p.cell].';
  E   = [state.p.E].';
  E_cell = accumarray(cid, E, [Nc 1], @sum, 0);

  if nargout>1
    stats = struct();
    stats.Nc = Nc;
    stats.Np = numel(state.p);
    stats.Nsp_cell = accumarray(cid, 1, [Nc 1], @sum, 0);
    stats.E_tot = sum(E_cell);
  end
end

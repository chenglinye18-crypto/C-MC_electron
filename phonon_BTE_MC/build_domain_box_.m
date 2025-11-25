function dom = build_domain_box_(mesh)
%BUILD_DOMAIN_BOX_  由 mesh 求计算域包围盒 [xmin xmax ymin ymax zmin zmax]
  if isstruct(mesh) && isfield(mesh,'domain_box') && numel(mesh.domain_box)==6
      dom = mesh.domain_box(:).';
      return;
  end
  if isstruct(mesh) && isfield(mesh,'boxes') && ~isempty(mesh.boxes)
      bx = mesh.boxes;
      dom = [min(bx(:,1)) max(bx(:,2)) ...
             min(bx(:,3)) max(bx(:,4)) ...
             min(bx(:,5)) max(bx(:,6))];
      return;
  end
  error('build_domain_box_: 缺少 mesh.domain_box 或 mesh.boxes，无法确定计算域。');
end

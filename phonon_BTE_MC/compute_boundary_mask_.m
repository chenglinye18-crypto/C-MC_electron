function mesh = compute_boundary_mask_(mesh, n_layers)
% not_boundary(c)=true 表示内部（PB禁用）；=false 表示靠外表面（PB启用）
  if nargin<2, n_layers = 1; end
  dom = build_domain_box_(mesh);
  bx  = mesh.boxes; Nc = size(bx,1); tol = 1e-12;

  touch = (abs(bx(:,1)-dom(1))<tol) | (abs(bx(:,2)-dom(2))<tol) | ...
          (abs(bx(:,3)-dom(3))<tol) | (abs(bx(:,4)-dom(4))<tol) | ...
          (abs(bx(:,5)-dom(5))<tol) | (abs(bx(:,6)-dom(6))<tol);
  bd = touch;                         % 触域者为边界层

  % 向里膨胀 (n_layers-1) 层
  for L = 2:n_layers
    grow = bd;
    for i = 1:Nc
      if bd(i), continue; end
      overlap = (bx(i,1)<=bx(:,2)+tol & bx(i,2)>=bx(:,1)-tol) & ...
                (bx(i,3)<=bx(:,4)+tol & bx(i,4)>=bx(:,3)-tol) & ...
                (bx(i,5)<=bx(:,6)+tol & bx(i,6)>=bx(:,5)-tol);
      if any(bd & overlap), grow(i)=true; end
    end
    bd = grow;
  end

  mesh.not_boundary = ~bd;   % 内部=1，边界层=0
end

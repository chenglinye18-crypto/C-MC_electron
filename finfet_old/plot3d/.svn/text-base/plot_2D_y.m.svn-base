function plot_2D_y(step,idx,ypos,x1,x2,z1,z2,ptype)
  function xend = find_index(x,a) 
        num = length(x);
        xend = num;
        for j = 1: num
             if (x(j) >= a) 
                xend = j;
                break;
             end
        end
  end
[x,y,z] = read_grid('../lgrid.txt');
numx = length(x);
numy = length(y);
numz = length(z);

if (ptype == 1) 
    name = strcat('../data/Electron',int2str(step));
else
    name = strcat('../data/Hole',int2str(step));
end

infile1=fopen(name);
flag = 1;
if (idx == 2) 
    scale = 1;
elseif (idx == 5)
    scale = 1e-6;
else 
    scale = 1;
end;
while flag
    id = fscanf(infile1, '%d', 3);
    if (numel(id) > 0)
       tmp = fscanf(infile1, '%f', 5);
       V(id(1) + 1, 1 + id(2), 1+id(3)) = tmp(idx); 
       if ((id(2) == 1)&&(tmp(4) > 0.1))
	 x(id(1) + 1)
	 z(id(3) + 1)
	 y(id(2) + 1)
	 tmp(5)
       end

    else 
       flag = 0;
    end
end

j0 = find_index(y, ypos);

for i = 1 : numx
    for k = 1 : numz
        data(i,k) = abs(V(i,j0,k)) * scale;
    end
end

xl(1) = x1;
xl(2) = x2;
yl(1) = y(1);
yl(2) = y(numy);
zl(1) = z1;
zl(2) = z2;

xbeg = find_index(x, xl(1));
xend = find_index(x, xl(2));
zbeg = find_index(z, zl(1));
zend = find_index(z, zl(2));

data(xbeg:xend,zbeg:zend);
surf(z(zbeg:zend),x(xbeg:xend),data(xbeg:xend,zbeg:zend));

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('z (\mu m)','FontSize', 30);
ylabel('x (nm)','FontSize', 30);
zlabel('Electron density (cm^{-3})', 'FontSize', 30);
%zlabel('Velocity (cm/s)', 'FontSize', 30);
%title('Electron density at z =25 nm', 'FontSize', 30);
%set(gca,'YDir','reverse')
view(125,25)
%set(gca,'XTick', [-0.085,-0.035,0,0.035,0.085]);
hold on
end

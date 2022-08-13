function plot_2D_x(step,idx,xpos,y1,y2,z1,z2,ptype)
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
    scale = 100;
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
    else 
       flag = 0;
    end
end

x0 = find_index(x, xpos);

for i = 1 : numy
    for k = 1 : numz
        data(i,k) = abs(V(x0,i,k)) * scale;
    end
end

xl(1) = x(1);
xl(2) = x(numx);
yl(1) = y1;
yl(2) = y2;
zl(1) = z1;
zl(2) = z2;

zbeg = find_index(z, zl(1));
zend = find_index(z, zl(2));
ybeg = find_index(y, yl(1));
yend = find_index(y, yl(2));

data(ybeg:yend,zbeg:zend);
surf(z(zbeg:zend),y(ybeg:yend),data(ybeg:yend,zbeg:zend));  

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('z (\mu m)','FontSize', 30);
ylabel('y (nm)','FontSize', 30);
zlabel('Electron density (cm^{-3})', 'FontSize', 30);
%zlabel('Velocity (cm/s)', 'FontSize', 30);
%title('Electron density at z =25 nm', 'FontSize', 30);
%set(gca,'YDir','reverse')
view(100,15)
%set(gca,'XTick', [-0.085,-0.035,0,0.035,0.085]);
hold on
end

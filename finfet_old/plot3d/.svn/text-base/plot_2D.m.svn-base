function plot_2D(step,idx,x1,x2,y1,y2,ptype)
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

[x,y] = read_grid('../grid.txt');
numx = length(x);
numy = length(y);

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
    id = fscanf(infile1, '%d', 2);
    if (numel(id) > 0)
       tmp = fscanf(infile1, '%f', 4);
       V(id(1) + 1, 1 + id(2)) = tmp(idx); 
    else 
       flag = 0;
    end
end

for i = 1 : numx
    for j = 1 : numy
        data(i,j) = abs(V(i,j)) * scale;
    end
end

xl(1) = x1;
xl(2) = x2;
yl(1) = y1;
yl(2) = y2;

xbeg = find_index(x, xl(1));
xend = find_index(x, xl(2));
ybeg = find_index(y, yl(1));
yend = find_index(y, yl(2));

data(xbeg:xend,ybeg:yend);
surf(y(ybeg:yend),x(xbeg:xend),data(xbeg:xend,ybeg:yend));  

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
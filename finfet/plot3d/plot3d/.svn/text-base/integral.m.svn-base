function integral(step,idx,y1,y2,x1,x2,z1,z2,ptype)
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

if (ptype == 1) 
    name = strcat('../qcdata/Electron',int2str(step));
else
    name = strcat('../qcdata/Hole',int2str(step));
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
       V(id(1) + 1, 1 + id(2), 1+id(3)) = tmp(idx) * scale; 
       if (idx == 5)
	 weight(id(1) + 1, id(2) + 1, id(3) + 1) = 1;
       else 
	 weight(id(1) + 1, id(2) + 1, id(3) + 1) = tmp(5);
       end
    else 
       flag = 0;
    end
end
flag = 1;
name = '../data/pvolume';
infile2=fopen(name);
while flag
    id = fscanf(infile2, '%d', 3);
    if (numel(id) > 0)
       tmp = fscanf(infile2, '%f', 1);
       area(id(1) + 1, id(2) + 1, id(3) + 1) = tmp;
    else 
       flag = 0;
    end
end

xl(1) = x1;
xl(2) = x2;
yl(1) = y1;
yl(2) = y2;
zl(1) = z1;
zl(2) = z2;

xbeg = find_index(x, xl(1));
xend = find_index(x, xl(2));
zbeg = find_index(z, zl(1));
zend = find_index(z, zl(2));
ybeg = find_index(y, yl(1));
yend = find_index(y, yl(2));

for j = ybeg : yend 
  data(j) = 0;
  tot_weight(j) = 0;
  for i = xbeg : xend
    for k = zbeg : zend
      tot_weight(j) = tot_weight(j) + area(i,j,k) * weight(i,j,k);
      data(j) = data(j) + V(i,j,k) * area(i,j,k) * weight(i,j,k);
    end
  end
  data(j) = data(j) / tot_weight(j);
end

plot(y(ybeg:yend),data(ybeg:yend),'k+-', 'LineWidth', 4,'MarkerSize',10);

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('y (nm)','FontSize', 30);
ylabel('Electron density (cm^{-3})','FontSize', 30);
%zlabel('Velocity (cm/s)', 'FontSize', 30);
%title('Electron density at z =25 nm', 'FontSize', 30);
%set(gca,'YDir','reverse')
%set(gca,'XTick', [-0.085,-0.035,0,0.035,0.085]);
hold on
end

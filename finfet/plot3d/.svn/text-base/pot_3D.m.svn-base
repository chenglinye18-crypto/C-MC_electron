function pot_3D(step)
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
name = strcat('../data/output/grid.txt');
infile1 = fopen(name);

numx = fscanf(infile1, '%d', 1);
for i = 1 : numx 
    x(i) = fscanf(infile1,'%f',1);
end
numy = fscanf(infile1,'%d', 1);
for i = 1 : numy
    y(i) = fscanf(infile1,'%f',1);
end
numz = fscanf(infile1,'%d', 1);
for i = 1 : numz
    z(i) = fscanf(infile1,'%f',1);
end
fclose(infile1);
name = strcat('../data/output/data/data/potential_',int2str(step));
infile1=fopen(name);
flag = 1;

while flag
    id = fscanf(infile1, '%d', 3);
    if (numel(id) > 0)
       tmp = fscanf(infile1, '%f', 1);
       V(id(2) + 1, 1 + id(3), 1+id(1)) = tmp; 
    else 
       flag = 0;
    end
end
x = x * 1000;
y = y * 1000;
z = z * 1000;

xl(1) = 1;
xl(2) = 161;
yl(1) = 0;
yl(2) = 40;
zl(1) = 0;
zl(2) = 40;

xbeg = find_index(x, xl(1));
xend = find_index(x, xl(2));
ybeg = find_index(y, yl(1));
yend = find_index(y, yl(2));
zbeg = find_index(z, zl(1));
zend = find_index(z, zl(2));

[xx,yy,zz] = meshgrid(z(zbeg:zend),y(ybeg:yend),x(xbeg:xend));

data = smooth3(V(ybeg:yend,zbeg:zend,xbeg:xend));
size(xx)
size(yy)
size(zz)
size(data)
h = slice(xx,yy,zz,data,[],[y(ybeg:yend)],[]);

box on
set(gca,'FontSize',20);
set(gca,'LineWidth',0.1);
%set(gca,'Box','on');
xlabel('z (\mu m)','FontSize', 30);
ylabel('y (\mu m)','FontSize', 30);
zlabel('x (\mu m)', 'FontSize', 30);
title('potential (V) ', 'FontSize', 30);
set(gca,'ZDir','reverse');
colorbar('EastOutside');
view(40,10);
%daspect([1 0.8 1]);
hold on
end
function plot_cdata_y(step,ypos,x1,x2,z1,z2,fid)

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

switch fid 
case 1
name = strcat('../data/xfield',int2str(step));
case 2
name = strcat('../data/yfield',int2str(step));
case 3
name = strcat('../data/zfield',int2str(step));
case 4
name = strcat('../data/electron_num',int2str(step));
end

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

xl(1) = x1;
xl(2) = x2;
zl(1) = z1;
zl(2) = z2;

zbeg = find_index(z, zl(1));
zend = find_index(z, zl(2));
xbeg = find_index(x, xl(1));
xend = find_index(x, xl(2));

i0 = find_index(y,ypos)

for i = 1 : numx - 1
    for k = 1 : numz - 1
        data(i,k) = V(i0,k,i);
    end
end

surf(z(zbeg:zend),x(xbeg:xend),data(xbeg:xend,zbeg:zend));

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('z / \mu m','FontSize', 30);
ylabel('x / \mu m','FontSize', 30);
zlabel('V', 'FontSize', 30);
%title('potential at z =25nm', 'FontSize', 30);
%set(gca,'YTick', 0:0.002:0.01);
hold on
end

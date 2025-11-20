function pot_2D_x(step,xpos,y1,y2,z1,z2, fid)

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
name = strcat('../data/qc',int2str(step));
case 2
name = strcat('../data/poipot',int2str(step));
case 3
name = strcat('../data/density',int2str(step));
case 4
name = strcat('../data/qpot',int2str(step));
case 5
name = strcat('../data/nq',int2str(step));
case 6
name = strcat('../data/ep',int2str(step));
case 7
name = strcat('../data/pot',int2str(step));
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

i0 = find_index(x,xpos);

for i = 1 : numy
    for k = 1 : numz
        data(i,k) = V(i,k,i0);
    end
end

surf(z(zbeg:zend),y(ybeg:yend),data(ybeg:yend,zbeg:zend));

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('z / \mu m','FontSize', 30);
ylabel('y / \mu m','FontSize', 30);
zlabel('V', 'FontSize', 30);
%title('potential at z =25nm', 'FontSize', 30);
%set(gca,'YTick', 0:0.002:0.01);
hold on
end

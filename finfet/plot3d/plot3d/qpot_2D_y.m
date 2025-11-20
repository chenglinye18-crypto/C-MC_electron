function qpot_2D_y(step,ypos,x1,x2,z1,z2)

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

[x,y,z] = read_grid('../grid.txt');

numx = length(x);
numy = length(y);
numz = length(z);

name1 = strcat('../data/QuanCorrection_',int2str(step));

infile1=fopen(name1);

name2 = strcat('../data/potential_',int2str(step));

infile2 = fopen(name2);

flag = 1;
while flag
    id = fscanf(infile1, '%d', 3);
    id1 = fscanf(infile2, '%d', 3);
    if (numel(id) > 0)
       tmp = fscanf(infile1, '%f', 1);
       tmp1 = fscanf(infile2, '%f', 1);
       V(id(2) + 1, 1 + id(3), 1+id(1)) = tmp + tmp1; 
    else 
       flag = 0;
    end
end

yl(1) = y(1);
yl(2) = y(numy);

xbeg = find_index(x, x1);
xend = find_index(x, x2);
zbeg = find_index(z, z1);
zend = find_index(z, z2);

j0 = find_index(y,ypos);

for i = 1 : numx
    for k = 1 : numz
        data(i,k) = V(j0,k,i) ;%* 0.025852;
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
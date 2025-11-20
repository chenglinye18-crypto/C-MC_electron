function pot_2D_z(step,zpos,x1,x2,y1,y2, fid)

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

xl(1) = x1;
xl(2) = x2;
yl(1) = y1;
yl(2) = y2;

xbeg = find_index(x, xl(1));
xend = find_index(x, xl(2));
ybeg = find_index(y, yl(1));
yend = find_index(y, yl(2));

k0 = find_index(z,zpos);

for i = 1 : numx
    for j = 1 : numy
        data(i,j) = V(j,k0,i);
    end
end

h = surf(y(ybeg:yend),x(xbeg:xend),data(xbeg:xend,ybeg:yend));
%set(h,'FaceColor', 'interp','EdgeColor','none');

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');

view(-40,25); 
text(-50,50,-1.3,'x (nm)', 'FontSize',40);
text(-5,157,-1.2,'y (nm)', 'FontSize',40);
text(33,-8,1.5, 'Potential along z = 5nm','FontSize',38);
text(38, 5, 1.0, 'V_g=1V,V_d=0.6V, with EP', 'FontSize', 38);

zlabel('Potential (V)', 'FontSize', 30);
set(gca,'YDir','reverse')
%title('potential at z =25nm', 'FontSize', 30);
%set(gca,'YTick', 0:0.002:0.01);
hold on
end

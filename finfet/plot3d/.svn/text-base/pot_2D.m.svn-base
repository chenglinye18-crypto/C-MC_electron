function pot_2D(step,x1,x2,y1,y2)
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

[x,y] = read_grid('../lgrid.txt');
numx = length(x);
numy = length(y);

%name = strcat('../data/QuanCorrection_',int2str(step));
%name = strcat('../data/qc',int2str(step));
%name = strcat('../data/potential_',int2str(step));
%name = strcat('../data/par.txt');
%name = strcat('../data/nq',int2str(step));
%name = strcat('../data/fermi_level',int2str(step));
%name = strcat('../data/pot',int2str(step));
name = strcat('../data/poipot',int2str(step));
%name = strcat('../data/base_poisson',int2str(step));
%name = strcat('../data/totpot',int2str(step));
%name = strcat('../data/density',int2str(step));
%name = strcat('../data/density-before',int2str(step));
%name = strcat('../data/density-after',int2str(step));
%name = strcat('../data/rhs',int2str(step));
infile1=fopen(name);
flag = 1;
while flag
    id = fscanf(infile1, '%d', 2);
    if (numel(id) > 0)
       tmp = fscanf(infile1, '%f', 1);
       V(id(2) + 1, 1+id(1)) = tmp; 
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

for i = 1 : numx    
    for k = 1 : numy
        data(i,k) = V(k,i); %* 0.025852;
    end
end

surf(y(ybeg:yend),x(xbeg:xend),data(xbeg:xend,ybeg:yend));

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('y / \mu m','FontSize', 30);
ylabel('x / \mu m','FontSize', 30);
%set(gca,'YTick',[0 10 20 30 40]);
zlabel('quantum correction (V)', 'FontSize', 30);
set(gca, 'YDir','reverse');

%title('potential at z =25nm', 'FontSize', 30);
%set(gca,'YTick', 0:0.002:0.01);
hold on
view(15,15);
end

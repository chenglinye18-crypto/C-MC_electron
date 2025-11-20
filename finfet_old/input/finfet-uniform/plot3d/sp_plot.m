function sp_plot(step,ypos,fileid)
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

if (fileid == 1)
name = strcat('../data/sp_dens',int2str(step));
else 
name = strcat('../data/sp_pot',int2str(step));
end

infile1=fopen(name);

numxz = fscanf(infile1, '%d', 2);

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

yid = find_index(y, ypos);

for i = 1 : numxz(1)    
    for k = 1 : numxz(2)
        data(i,k) = V(i,k,yid);
    end
end

surf(1:numxz(1), 1:numxz(2), data);

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

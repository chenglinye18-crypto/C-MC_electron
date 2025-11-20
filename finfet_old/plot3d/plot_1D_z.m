function plot_1D_z(step,idx,z1,z2, xpos,ypos, ptype)
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

if (ptype == 1) 
    name = strcat('../data/Electron_',int2str(step));
else
    name = strcat('../data/Hole_',int2str(step));
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
    
i0 = find_index(x, xpos);
j0 = find_index(y, ypos);

for k = 1 : numz
  data(k) = abs(V(i0,j0,k)) * scale;
end

zbeg = find_index(z, z1);
zend = find_index(z, z2);

plot(z(zbeg:zend), data(zbeg:zend),'r');

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('z (\mu m)','FontSize', 30);
%ylabel('y (nm)','FontSize', 30);
ylabel('Electron density (cm^{-3})', 'FontSize', 30);
%zlabel('Velocity (cm/s)', 'FontSize', 30);
%title('Electron density at z =25 nm', 'FontSize', 30);
%set(gca,'YDir','reverse')
%view(105,25)
%set(gca,'XTick', [-0.085,-0.035,0,0.035,0.085]);
hold on
end

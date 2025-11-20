function pot_1D_z(step,z1,z2,ipos,jpos)

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

%name = strcat('../data/QuanCorrection_',int2str(step));
name = strcat('../data/potential_',int2str(step));
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

zbeg = find_index(z, z1)
zend = find_index(z, z2)

j0 = find_index(y,jpos)
i0 = find_index(x,ipos)

for i = 1 : numz    
  data(i) = V(j0,i,i0) ;%* 0.025852;
end
data
plot(z(zbeg:zend), data(zbeg:zend));
set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('z / \mu m','FontSize', 30);
%ylabel('x / \mu m','FontSize', 30);
ylabel('V', 'FontSize', 30);
%title('potential at z =25nm', 'FontSize', 30);
%set(gca,'YTick', 0:0.002:0.01);
hold on
end

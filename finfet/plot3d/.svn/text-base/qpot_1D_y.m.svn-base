function qpot_1D_y(step,y1,y2,ipos,kpos)

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

name = strcat('../data/QuanCorrection_',int2str(step));
infile1=fopen(name);
name2 = strcat('../data/potential_',int2str(step));
infile2=fopen(name2);
flag = 1;

while flag
    id = fscanf(infile1, '%d', 3);
    id1 = fscanf(infile2, '%d', 3);
    if (numel(id) > 0)
       tmp = fscanf(infile1, '%f', 1);
       tmp1 = fscanf(infile2, '%f',1);
       V(id(2) + 1, 1 + id(3), 1+id(1)) = tmp + tmp1; 
    else 
       flag = 0;
    end
end

ybeg = find_index(y, y1);
yend = find_index(y, y2);

i0 = find_index(x,ipos)
k0 = find_index(z,kpos)

for i = 1 : numy    
  data(i) = V(i,k0,i0) ;%* 0.025852;
end

plot(y(ybeg:yend), data(ybeg:yend),'m');

set(gca,'FontSize',25);
set(gca,'LineWidth',3);
%set(gca,'Box','on');
xlabel('y / \mu m','FontSize', 30);
%ylabel('x / \mu m','FontSize', 30);
ylabel('V', 'FontSize', 30);
%title('potential at z =25nm', 'FontSize', 30);
%set(gca,'YTick', 0:0.002:0.01);
hold on
end

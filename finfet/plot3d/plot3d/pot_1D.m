function pot_1D(step,ix,k0)
name = strcat('/home/zhangwei/study/tmp/hpc_data/3D/output/grid.txt');
%name = strcat('/home/zhangwei/study/tmp/hpc_data/3D/FinFet/grid.txt');
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
%x
%y
%name = strcat('/home/zhangwei/study/tmp/hpc_data/history_data/diff_voltage/0.2/output/statistic/Electron/', int2str(step));
%name = strcat('/home/zhangwei/study/tmp/hpc_data/2D/output/statistic/Hole/', int2str(step));
%name = strcat('/home/zhangwei/study/tmp/hpc_data/tmp/data/output/Electron_', int2str(step));
name = strcat('/home/zhangwei/study/tmp/hpc_data/3D/output/data/potential_',int2str(step));
%name = strcat('/home/zhangwei/study/tmp/hpc_data/3D/FinFet/Potential',int2str(step));
%name = strcat('/home/zhangwei/study/tmp/hpc_data/tmp/data/output1/Hole_', int2str(step));
infile1=fopen(name);
flag = 1;

while flag
    id = fscanf(infile1, '%d', 3);
    if (numel(id) > 0)
       tmp = fscanf(infile1, '%f', 1);
       V(id(1) + 1, 1 + id(2), 1+id(3)) = tmp;% *  0.025852; 
    else 
       flag = 0;
    end
end
for k = 1 : numy
   data(k) = V(ix,k,k0);
end

plot(y,data,'r');

hold on


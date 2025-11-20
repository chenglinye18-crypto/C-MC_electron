function [x,y,z] = read_grid(gname)

infile1 = fopen(gname);

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

x = x * 1000;
y = y * 1000;
z = z * 1000;

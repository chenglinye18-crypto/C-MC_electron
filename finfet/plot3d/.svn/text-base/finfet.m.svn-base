function finfet
name = strcat('/home/zhangwei/study/tmp/hpc_data/3D/FinFet/grid.txt');
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
% 
vert = [0 1 1;0 1 11;0 1 41;75 1 1;75 1 41;75 11 41;75 11 1; 0 11 1;75 1 11;75 11 11];
%vert1 = []; 
fac = [1 2 9 4;2 3 5 9;9 5 6 10;4 9 10 7;1 4 7 8];
tcolor = [1 0 0;0.7 0.7 0.7;0.7 0.7 0.7;1 0 0;1 0 0];
view(10,10);
set(gca,'ZDir','reverse');
patch('Faces' , fac, 'Vertices', vert, 'FaceVertexCData', tcolor, 'FaceColor', 'flat');

vert1 = [25 0 0;25 0 12;50 0 12;50 0 0;25 12 0; 50 12 0];
fac1 = [1 2 3 4; 1 4 6 5];
tcolor1 = [0 1 0.8];
patch('Faces' , fac1, 'Vertices', vert1, 'FaceVertexCData', tcolor1, 'FaceColor', 'flat');

vert2 = [50 0 12;50 1 12;50 1 1; 50 11 1; 50 11 12;50 12 12; 50 12 0;50 0 0;];
fac2 = [1 2 3 4 5 6 7 8];
patch('Faces' , fac2, 'Vertices', vert2, 'FaceVertexCData', tcolor1, 'FaceColor', 'flat');

axis tight
set(gca,'FontSize',25);
set(gca,'LineWidth',3);
text(65, 0, 0,'Drain','FontSize',30,'Color','b');
text(10, 0, 0, 'Source','FontSize',30,'Color','b');
text(37, 0, 0, 'Gate','FontSize',30,'Color','b');
text(37, 0, 0 , 'Substrate','FontSize',30,'Color','b');

xlabel('y (nm)','FontSize', 30);
ylabel('z (nm)','FontSize', 30);
zlabel('x (nm)','FontSize', 30);
function trace(par_id)

name = strcat('../data/trace',int2str(par_id));
infile1=fopen(name);

flag = 1;
while flag
    id = fscanf(infile1, '%d', 1);
    if (numel(id) > 0)
       cell(id + 1,:) = fscanf(infile1, '%d', 3);
       pos(id + 1,:) = fscanf(infile1, '%f', 4);
       pos(id + 1,:) = pos(id + 1,:) * 1e9;
    else 
       flag = 0;
    end
end


y=[-30 30];
x = [0 10];
z = [0 10];


line([0 0], y, [0 0]);
line([10 10], y, [0 0]);
line([0 0], y, [10 10]);
line([10 10], y, [10 10]);
line([0 10], [-30 -30], [10 10]);
line([0 10], [-30 -30], [0 0]);
line([0 10], [30 30], [10 10]);
line([0 10], [30 30], [0 0]);

line([0 10], [-10 -10], [10 10]);
line([0 10], [-10 -10], [0 0]);
line([0 10], [10 10], [10 10]);
line([0 10], [10 10], [0 0]);

line([0 0], [-30 -30], [0 10]);
line([10 10], [-30 -30], [0 10]);
line([0 0], [30 30], [0 10]);
line([10 10], [30 30], [0 10]);

line([0 0], [-10 -10], [0 10]);
line([10 10], [-10 -10], [0 10]);
line([0 0], [10 10], [0 10]);
line([10 10], [10 10], [0 10]);

view(90,0)
%comet3(pos(:,1), pos(:,1), pos(:,1));


ylabel('y', 'FontSize', 30);
hold on

comet3(pos(:,1), pos(:,2), pos(:,3));

pos(:,4)


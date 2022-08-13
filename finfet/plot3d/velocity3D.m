function velocity3D(step,idx,x1,x2,z1,z2)
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

name = strcat('../data/Electron_',int2str(step));
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
       V(id(2) + 1, 1 + id(3), 1+id(1)) = abs(tmp(idx)) * scale; 
    else 
       flag = 0;
    end
end

xl(1) = x1;
xl(2) = x2;
zl(1) = z1;
zl(2) = z2;


xbeg = find_index(x, xl(1));
xend = find_index(x, xl(2));
zbeg = find_index(z, zl(1));
zend = find_index(z, zl(2));
ybeg = find_index(y, 50);
yend = find_index(y, 70);

[xx,yy,zz] = meshgrid(z(zbeg:zend),y(ybeg:yend),x(xbeg:xend));

data = smooth3(V(ybeg:yend, zbeg:zend,xbeg:xend));

h = slice(xx,yy,zz,data,[],[50:1:55],[]);

sh = get(h);
%set(h(1),'FaceAlpha',0.3);
%set(h(1),'AlphaDataMapping','none');
%set(h(1),'AlphaData',0.2);
for i = 1 : 0 %length(h)
set(h(i),'EdgeAlpha',0);
set(h(i),'EdgeColor',[0 0 0]);
set(h(i),'FaceColor','flat');
set(h(i),'LineStyle','none');
end
%sh = get(h)
%sh(1).AlphaData
%sh(1).AlphaDataMapping
%sh(1).CData
%set(h(1),'FaceColor','none');

%p = patch(get(h));
%set(p,'facecolor','red','edgecolor','none');
%camlight; lighting gouraud;
%alpha(.5);
%alpha('color')
%set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
%alphamap('rampdown')
%alphamap('increase',.2)
%alpha(0.5);
%colormap(hsv)
%s=volumeVisualization(xx,yy,zz, V(:,4:23,xbeg:xend));
%s.addSlicePlane(0.035)
set(gca,'ZDir','reverse');
%et(gca,'XDir','reverse')
view(40,10); daspect([1 0.25 1]);
axis tight
%view(110,12); 
%daspect([1 1 1]);axis tight
set(gca,'FontSize',40);
set(gca,'LineWidth',5);
%text(z(2) + 0.01 ,y(74), -0.001 ,'Drain','FontSize',30,'Color','y');
%text(z(2) +0.01,y(6),  -0.001, 'Source','FontSize',30,'Color','y');

%line([z(3) z(3)], [0.035 0.035],[x(xbeg),x(xend)],'Color','k','LineWidth',2);
%line([z(1) z(3)], [0.035 0.035],[x(xbeg),x(xbeg)],'Color','k','LineWidth',2);
%line([z(3) z(3)], [-0.035 -0.035],[x(xbeg),x(xend)],'Color','k','LineWidth',2);
%line([z(1) z(3)], [-0.035 -0.035],[x(xbeg),x(xbeg)],'Color','k','LineWidth',2);
%set(gca,'Box','on');
xlabel('z (nm)','FontSize', 40);
ylabel('y (nm)','FontSize', 40);
zlabel('x (nm)','FontSize', 40);
%title('Electron velocity (cm/s)', 'FontSize', 45); 
title('Electron density (cm^{-3})', 'FontSize', 45);
%title('', 'FontSize', 30);
%title('Electron density at z =25 nm', 'FontSize', 30);
%set(gca,'YDir','reverse')
%view(110,10)
%set(gca,'ZTick', [1,5,11]);
set(gca,'XTick',[15:5:25]);
set(gca,'ZTick',[1 :2: 10]);
set(gca,'YTick',[50:1:55]);
%set(gca,'XTick', [0,0.05]);
%set(gca,'TickLength',[0.01,0.01]);
%colorbar('horiz');
h = colorbar('EastOutside')
%set(get(h,'Title'),'String','cm/s','FontSize', 35);
set(h,'FontSize',30,'FontWeight', 'bold');
map = colormap;
N = length(map);
maxv = 1e20;
minv = 1e19;
for i = 1 : N
    linear_dis = minv + (maxv - minv) * (i - 1) / (N - 1);
    log_index = round((log(linear_dis) - log (minv)) / (log(maxv) - log(minv)) * (N - 1))
    for k = 1 : 3
        newmap(i,k) = map(1 + log_index, k);
    end
end
colormap(newmap);
%pbaspect([1 10 2])
%rotate(s.h,[1,0,0],45);
hold on
end
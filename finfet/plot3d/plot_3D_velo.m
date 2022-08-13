function plot_3D_velo(step,idx,x1,x2,y1,y2,z1,z2,slice_num)

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

name = strcat('../data/Electron',int2str(step));

infile1=fopen(name);
flag = 1;
if (idx == 2) 
    scale = 100;
elseif (idx == 5)
    scale = 1e-26;
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

xbeg = find_index(x, x1);
xend = find_index(x, x2);
ybeg = find_index(y, y1);
yend = find_index(y, y2);
zbeg = find_index(z, z1);
zend = find_index(z, z2);

[xx,yy,zz] = meshgrid(z(zbeg:zend),y(ybeg:yend),x(xbeg:xend));
%z(24)
%x(14)
%size(xx)
%size(yy)
%size(zz)
%slice(xx,yy,zz,V(:,4:25,xbeg:xend),[z(4),z(24)], 0.025:0.005:0.05,x(5));
     
data = smooth3(V(ybeg:yend,zbeg:zend,xbeg:xend));

max(max(max(data)))

%posy = [-9,-7,-5,-3,-1];
%posy = [1,3,5,7,9];
posy = [-9,-7,-5,-3, -1,1,3,5,7,9];

hslices = slice(xx,yy,zz,data,[],posy,[],'nearest');
sh = get(hslices);


%for i = 1 : slice_num
%set(hslices(i),'FaceAlpha',0.8);
%set(hslices(1),'AlphaDataMapping','none');
%set(hslices(1),'AlphaData',0.2);
%set(hslices(i),'EdgeAlpha',0.2);
%set(hslices(i),'EdgeColor',[0 0 0]);
%set(hslices(i),'FaceColor','flat');
%set(hslices(i),'LineStyle','none');
%end

set(hslices,'FaceColor', 'interp','EdgeColor','none');
camzoom(1.4);
%box off
grid off

set(gca,'ZDir','reverse');
%set(gca,'XDir','reverse');
pos = get(gca,'Position');
pos(3) = pos(3) * 0.85;
set(gca,'Position',pos);
view(50,5); daspect([1 0.5 3]);
%axis tight
set(gca,'FontSize',35);
set(gca,'LineWidth',3);
xlim([-1,11]);
zlim([-1,41]);
ylim([-10,9.5]);
set(gca,'ZTick',[0,5,10,15,20,25,30,35,40]);
set(gca,'XTick',[0,2,4,6,8,10]);
set(gca,'YTick',[-9,-7,-5,-3,-1,1,3,5,7,9]);
%set(gca,'YTick',[-9,-7,-5,-3,-1,1,3,5,7,9]);
%set(gca,'YTick',[-8,-6,-4,-2,0,2,4,6,8]);
%set(gca,'YTick',[-8,-6,-4,-2,0,2,4,6,8]);
set(gca,'TickDir','in');
set(gca,'TickLength',[0.015,0.015]);
%xlabel('z (nm)','FontSize', 40);
text(8,-13,50,'z (nm)', 'FontSize',40);
text(12,-1,50,'y (nm)', 'FontSize',40);
text(0,-9,-13.5, 'Velocity slices with EP (cm/s)', 'FontSize', 38);
text(0,-9,-3.5, 'V_g=1V,V_d=0.6V','FontSize',38);
%text(0,-8,-0.6,'V_g=V_d=1V','FontSize',38);
%ylabel('y (nm)','FontSize', 40);
zlabel('x (nm)','FontSize', 40);
%title('Electron density (cm^{-3})', 'FontSize', 50);
h= colorbar('EastOutside');
axes(h);
%color = colormap;

color = colormap(jet(64));

h3 = 80;

newcolor = color;

for i = 64 : h3
  newcolor(i,:) = color(64,:);
end

%for i = h1 : 64
%  j = round(h2 + (i - h1) / (64 - h1) * (64 - h2));
%  newcolor(i,:) = color(j,:);
%end

colormap(newcolor);


%set(get(h,'Title'),'String','10^{19} cm^{-3}','FontSize', 40);
get(h,'Title')
%caxis([10,24]);
%colormap((cool(24)));
set(h,'FontSize',40,'FontWeight', 'bold');
barpos = get(h,'Position');
barpos(1) = pos(1) + pos(3) * 1.075;
barpos(4) = barpos(4) * 0.95;
set(h,'Position',barpos);
%text(3, 1.2, 'x 10^{20}','FontSize',40);
%get(h,'XTickLabel',[2,4,6,8,10,12])
%set(h,'YLim',[0.01,0.7])
set(h,'YLim',[0.4e7,2.1e7]);
%set(h,'YTick',[0.3;0.6;0.9;1.2;1.5]);
%set(h,'YTick',[0.1;0.2;0.3;0.4;0.5;0.6]);
%set(h,'YTickLabel',[0.1;0.2;0.3;0.4;0.5;0.6]);
set(h,'YTick',[0.4;0.8;1.2;1.6;2.0] * 1e7);
set(h,'YTickLabel',[0.4;0.8;1.2;1.6;2.0]);
%text(barpos(1) + barpos(3) , barpos(4) * 0.9,'x 10^{19}','FontSize', 40);
hold on
end

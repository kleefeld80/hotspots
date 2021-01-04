% TODO
% normalize the eigenfunction caxis([0 1])
function plotthething5holes(filename,xx,yy,Z,outer,inner,inner2,inner3,inner4,inner5,fignr,xticksvec,yticksvec)
    terrain = [0 166 0
        22 174 0
        45 182 0
        71 190 0
        99 198 0
        128 206 0
        160 214 0
        194 222 0
        230 230 0
        231 208 26
        232 192 52
        234 182 78
        235 177 105
        237 178 131
        238 185 159
        239 198 186
        241 217 214
        242 242 242];
    terrain = terrain/255;
    
    figure(fignr)
       hold on
       ticks = 12; % Schriftgroesse der Achsenunterteilungen
       set(gca,'fontsize',ticks,'fontweight','bold','LineWidth',2);
        %Z = Z/Z(min(find(abs(Z(:)) == max(abs(Z(:))))));
        [C,h]=contourf(xx,yy,Z,40);
        %h.LineStyle ='none';
        line(outer(1,:),outer(2,:),'Color','k','Linewidth',4)
        line(inner(1,:),inner(2,:),'Color','k','Linewidth',4)
        line(inner2(1,:),inner2(2,:),'Color','k','Linewidth',4)
        line(inner3(1,:),inner3(2,:),'Color','k','Linewidth',4)
        line(inner4(1,:),inner4(2,:),'Color','k','Linewidth',4)
        line(inner5(1,:),inner5(2,:),'Color','k','Linewidth',4)
        cb=colorbar;
        cb.Ticks = linspace(0,1,5);
        colormap(terrain)
        c=colorbar;
        c.LineWidth=2;
        c.Color='black';
        [row1, col1] = find(ismember(Z, min(Z(:))));
        [row2, col2] = find(ismember(Z, max(Z(:))));
        xmin=xx(row1,col1);
        ymin=yy(row1,col1);
        xmax=xx(row2,col2);
        ymax=yy(row2,col2);
        plot(xmin,ymin,'ro','MarkerFaceColor','r')
        plot(xmax,ymax,'bo','MarkerFaceColor','b')
        xticks(xticksvec)
        yticks(yticksvec)
        size = 12; % Schriftgroesse der Achsenbeschriftung
        xlabel('x','FontSize',size,'fontweight','bold');
        ylabel('y','FontSize',size,'fontweight','bold');
        axis equal
        box on
        hold off
        print(filename,'-depsc2')
end
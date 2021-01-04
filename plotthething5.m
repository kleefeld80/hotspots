% TODO
% normalize the eigenfunction caxis([0 1])
function plotthething5(filename,outer,inner,fignr,xticksvec,yticksvec)
    
    figure(fignr)
       hold on
       ticks = 12; % Schriftgroesse der Achsenunterteilungen
       set(gca,'fontsize',ticks,'fontweight','bold','LineWidth',2);
        %Z = Z/Z(min(find(abs(Z(:)) == max(abs(Z(:))))));
        
        %h.LineStyle ='none';
        fill(outer(1,:),outer(2,:),[0.8,0.8,0.8],'Linewidth',4)
        fill(inner(1,:),inner(2,:),[1,1,1],'Linewidth',4)
        
        xticks(xticksvec)
        yticks(yticksvec)
        size = 12; % Schriftgroesse der Achsenbeschriftung
        xlabel('x','FontSize',size,'fontweight','bold');
        ylabel('y','FontSize',size,'fontweight','bold');
        axis equal
        axis([-1.8,1.8,-1.8,1.8])
        box on
        hold off
        print(filename,'-depsc2')
end
function plotTJ(MT, MD)

%% Trajetória Mínima
% Linha e Coluna da menor distância
linhaMin = min(MD);
valorMin = min(linhaMin);
[LMin, CMin] = find(MD==valorMin);
[Lmin,Cmin] = size(MT{LMin,CMin});

figure(15)
hold on
title('Trajetória Mínima')
for i=Lmin:Cmin
    if(i<Cmin)
%         title(['t=' num2str(i)]);
        plot(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        x = [MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i+1}(1)];
        y = [MT{LMin,CMin}{1,i}(2),MT{LMin,CMin}{1,i+1}(2)];
        line(x,y,'Color','red','LineWidth',1);
%         drawnow
    else
        plot(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        text(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),strcat('(',num2str(MT{LMin,CMin}{1,i}(1),2),',',num2str(MT{LMin,CMin}{1,i}(2),2),')'), ...
                'horiz','center', ...
                'vert','bottom', ...
                'FontWeight','bold', ...
                'EdgeColor','black')
%         drawnow
    end
end
% Texto destacado da primeira coordenada mínima
text(MT{LMin,CMin}{1,1}(1),MT{LMin,CMin}{1,1}(2),strcat('(',num2str(MT{LMin,CMin}{1,1}(1),2),',',num2str(MT{LMin,CMin}{1,1}(2),2),')'), ...
    'horiz','center', ...
    'vert','bottom', ...
    'FontWeight','bold', ...
    'EdgeColor','black')
axis([-1 1 -1 1])
hold off

%% Trajetória Máxima
% Linha e Coluna da maior distância
linhaMax = max(MD);
valorMax = max(linhaMax);
[LMax, CMax] = find(MD==valorMax);
[Lmax,Cmax] = size(MT{LMax,CMax});

figure(16)
hold on
title('Trajetória Máxima')
for i=Lmax:Cmax
    if(i<Cmax)
%         title(['t=' num2str(i)]);
        plot(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        x = [MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i+1}(1)];
        y = [MT{LMax,CMax}{1,i}(2),MT{LMax,CMax}{1,i+1}(2)];
        line(x,y,'Color','red','LineWidth',1);
%         drawnow
    else
        plot(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        text(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),strcat('(',num2str(MT{LMax,CMax}{1,i}(1),2),',',num2str(MT{LMax,CMax}{1,i}(2),2),')'), ...
                'horiz','center', ...
                'vert','bottom', ...
                'FontWeight','bold', ...
                'EdgeColor','black')
%         drawnow
    end
end
% Texto destacado da primeira coordenada máxima
text(MT{LMax,CMax}{1,1}(1),MT{LMax,CMax}{1,1}(2),strcat('(',num2str(MT{LMax,CMax}{1,1}(1),2),',',num2str(MT{LMax,CMax}{1,1}(2),2),')'), ...
    'horiz','center', ...
    'vert','bottom', ...
    'FontWeight','bold', ...
    'EdgeColor','black')
axis([-1 1 -1 1])
hold off

end

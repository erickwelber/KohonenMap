function plotTJ(MT, MD, x1, x2, w1, w2, varargin)

figure(8)

% configurações do plot
if strcmp(varargin, 'annotation')
    set(gca, 'FontName', 'Helvetica', 'FontSize', 12)
    axis([-1 1 -1 1])
    title('Trajetória Máxima e Mínima')
end

%% Mapa Auto-Organizado
% figure(8)
plot(x1,x2,'ob')
hold on
plot(w1,w2,'r','linewidth',2)
plot(w1',w2','r','linewidth',2)
plot(w1,w2,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
hold on
% title('Trajetória Máxima e Mínima')
%% Trajetória de todos os neurônios
% tamanho da matriz
% [M,N] = size(MT);
% figure(7)
% for m=1:M
%     for n=1:N
%         [L,C] = size(MT{m,n});
%         for i=L:C
%             if(i+1<=C)
%                 plot(MT{1,1}{1,i}(1),MT{1,1}{1,i}(2),'bo','LineWidth',2);
%                 x = [MT{1,1}{1,i}(1),MT{1,1}{1,i+1}(1)];
%                 y = [MT{1,1}{1,i}(2),MT{1,1}{1,i+1}(2)];
%                 line(x,y,'LineStyle','--','LineWidth',2);
%             else
%                 plot(MT{1,1}{1,i}(1),MT{1,1}{1,i}(2),'bo','LineWidth',2);
%             end
%             hold on
%         end
%     end
% end

%% Trajetória Mínima
% Linha e Coluna da menor distância
linhaMin = min(MD);
valorMin = min(linhaMin);
[LMin, CMin] = find(MD==valorMin);
[Lmin,Cmin] = size(MT{LMin,CMin});
% figure(8)
% title('Trajetória Máxima e Mínima')
for i=Lmin:Cmin
    if(i+1<=Cmin)
%         title(['t=' num2str(i)]);
        plot(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        x = [MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i+1}(1)];
        y = [MT{LMin,CMin}{1,i}(2),MT{LMin,CMin}{1,i+1}(2)];
        line(x,y,'Color','red','LineWidth',0.5);
%         drawnow
    else
        plot(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
%         drawnow
    end
    if(i==Cmin)
        text(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),strcat('(',num2str(MT{LMin,CMin}{1,i}(1),2),',',num2str(MT{LMin,CMin}{1,i}(2),2),')'), ...
                'horiz','center', ...
                'vert','bottom', ...
                'FontWeight','bold', ...
                'EdgeColor','black')
    end
    hold on
end

%% Trajetória Máxima
% Linha e Coluna da maior distância
linhaMax = max(MD);
valorMax = max(linhaMax);
[LMax, CMax] = find(MD==valorMax);
[Lmax,Cmax] = size(MT{LMax,CMax});
% figure(8)
% title('Trajetória Máxima e Mínima')
for i=Lmax:Cmax
    if(i+1<=Cmax)
%         title(['t=' num2str(i)]);
        plot(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        x = [MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i+1}(1)];
        y = [MT{LMax,CMax}{1,i}(2),MT{LMax,CMax}{1,i+1}(2)];
        line(x,y,'Color','red','LineWidth',0.5);
%         drawnow
    else
        plot(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
%         drawnow
    end
    if(i==Cmax)
        text(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),strcat('(',num2str(MT{LMax,CMax}{1,i}(1),2),',',num2str(MT{LMax,CMax}{1,i}(2),2),')'), ...
                'horiz','center', ...
                'vert','bottom', ...
                'FontWeight','bold', ...
                'EdgeColor','black')
    end
    hold on
end

end

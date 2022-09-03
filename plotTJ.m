function plotTJ(MT, MD, varargin)

% configurações do plot
if strcmp(varargin, 'annotation')
    set(gca, 'FontName', 'Helvetica', 'FontSize', 12)
    xlim([-1 1]);
    ylim([-1 1]);
end

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
figure(8)
for i=Lmin:Cmin
    if(i+1<=Cmin)
        title(['Trajetória Mínima=' num2str(i)]);
        plot(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        x = [MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i+1}(1)];
        y = [MT{LMin,CMin}{1,i}(2),MT{LMin,CMin}{1,i+1}(2)];
        line(x,y,'Color','red','LineWidth',2);
        drawnow
    else
        plot(MT{LMin,CMin}{1,i}(1),MT{LMin,CMin}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        drawnow
    end
    hold on
end

%% Trajetória Máxima
% Linha e Coluna da maior distância
linhaMax = max(MD);
valorMax = max(linhaMax);
[LMax, CMax] = find(MD==valorMax);
[Lmax,Cmax] = size(MT{LMax,CMax});
figure(9)
for i=Lmax:Cmax
    if(i+1<=Cmax)
        title(['Trajetória Mínima=' num2str(i)]);
        plot(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        x = [MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i+1}(1)];
        y = [MT{LMax,CMax}{1,i}(2),MT{LMax,CMax}{1,i+1}(2)];
        line(x,y,'Color','red','LineWidth',2);
        drawnow
    else
        plot(MT{LMax,CMax}{1,i}(1),MT{LMax,CMax}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
        drawnow
    end
    hold on
end

end

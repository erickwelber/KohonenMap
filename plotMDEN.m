function plotMDEN(MDEN)

% Plot de sinal do neurônio 1 com o neurônio 2 ao longo de 300 iterações
% [m,n] = size(MDEN{1,1});
% figure(10)
% hold on
% title('Visualização dos Neurônios 1 e 2')
% for i=m:n
%     if(i<n)
%         %Neurônio 1
%         plot(MDEN{1,1}{1,i}(1),MDEN{1,1}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
%         x1 = [MDEN{1,1}{1,i}(1),MDEN{1,1}{1,i+1}(1)];
%         y1 = [MDEN{1,1}{1,i}(2),MDEN{1,1}{1,i+1}(2)];
%         line(x1,y1,'Color','red','LineWidth',1);
%         
%         %Neurônio 2
%         plot(MDEN{1,1}{1,i}(3),MDEN{1,1}{1,i}(4),'bo','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','c');
%         x2 = [MDEN{1,1}{1,i}(3),MDEN{1,1}{1,i+1}(3)];
%         y2 = [MDEN{1,1}{1,i}(4),MDEN{1,1}{1,i+1}(4)];
%         line(x2,y2,'Color','blue','LineWidth',1);
%         
%     else
%         %Neurônio 1
%         plot(MDEN{1,1}{1,i}(1),MDEN{1,1}{1,i}(2),'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
%         
%         %Neurônio 2
%         plot(MDEN{1,1}{1,i}(3),MDEN{1,1}{1,i}(4),'bo','MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor','c');
%         
%     end
% end
% % linha de distância entre os neurônios 1 e 2 da primeira coordenada
%     xx1=[MDEN{1,1}{1,1}(1),MDEN{1,1}{1,1}(3)];
%     yy1=[MDEN{1,1}{1,1}(2),MDEN{1,1}{1,1}(4)];
%     plot(xx1,yy1,'Color','black','LineWidth',1);
%     text(xx1,yy1,strcat('[',num2str(MDEN{1,1}{1,1}(5),4),']'));
% axis([-1 1 -1 1])
% hold off
% legend('Neurônio 1','','Neurônio 2','')

% Plot de distâncias entre os neurônios
figure(11)
hold on
title('Distâncias Sinápticas')
[M,N] = size(MDEN);
legendInfo = {};
for i=1:M
    aux = 0;
    for j=1:N
        if(~isempty(MDEN{i,j}))
            data = [];
            TAM = length(MDEN{i,j});
            for ii=1:TAM
                data(ii) = MDEN{i,j}{1,ii}(5);
            end
            plot(300/length(data).*(1:1:length(data)),data)
            grid on
            aux = aux + 1;
            legendInfo{end+1} = ['Entre ' num2str(i) ' e ' num2str(i+aux)];
        end
    end
end
hold off
legend(legendInfo,'Location','northwest','NumColumns',2)
end
function plotCS(MD, MC, w1, w2, varargin)

% tamanho da matriz
[M,N] = size(MD);

% valor mínimo e índice
linhaMin = min(MD);
valorMin = min(linhaMin);
disp(valorMin);
[LMDMin, CMDMin] = find(MD==valorMin); % Linha e Coluna do menor valor na Matriz de Distância

% coordenada mínima de MC
vetorMin = cell2mat(MC(LMDMin,CMDMin)); % vetor com as respectivas coordenada do valor mínimo
disp(vetorMin);
valorCW1Min = vetorMin(1); % valor da coordenada do menor valor de w1
valorCW2Min = vetorMin(2); % valor da coordenada do menor valor de w2

% valor máximo e índice de MD
linhaMax = max(MD);
valorMax = max(linhaMax);
disp(valorMax);
[LMDMax, CMDMax] = find(MD==valorMax); % Linha e Coluna do maior valor na Matriz de Distância

% coordenada máxima de MC
vetorMax = cell2mat(MC(LMDMax,CMDMax));
disp(vetorMax);
valorCW1Max = vetorMax(1); % valor da coordenada do maior valor de w1
valorCW2Max = vetorMax(2); % valor da coordenada do maior valor de w2

% configurações de visualização da matriz
figure(6)
plot(w1,w2,'r','linewidth',2)
hold on
plot(w1',w2','r','linewidth',2)
hold on

% configurações do plot
if strcmp(varargin, 'annotation')
    % x-label, y-label, x-ticks, y-ticks, title
    set(gca, 'FontName', 'Helvetica', 'FontSize', 12)
%     xlabel('w1')
%     ylabel('w2')
    title('Coordenadas sinápticas')
end

for m = 1:M
    for n = 1:N
        if((w1(m,n)==valorCW1Min) && (w2(m,n)==valorCW2Min))
            figure(6)
            plot(w1(m,n),w2(m,n),'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            text(w1(m,n),w2(m,n),strcat('(',num2str(w1(m,n),2),',',num2str(w2(m,n),2),')'), ...
                'horiz','center', ...
                'vert','bottom', ...
                'FontWeight','bold', ...
                'EdgeColor','black')
            hold on
        elseif((w1(m,n)==valorCW1Max) && (w2(m,n)==valorCW2Max))
            figure(6)
            plot(w1(m,n),w2(m,n),'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            text(w1(m,n),w2(m,n),strcat('(',num2str(w1(m,n),2),',',num2str(w2(m,n),2),')'), ...
                'horiz','center', ...
                'vert','bottom',...
                'FontWeight','bold', ...
                'EdgeColor','black');
            hold on
        else
            figure(6)
            plot(w1(m,n),w2(m,n),'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
            text(w1(m,n),w2(m,n),strcat('(',num2str(w1(m,n),2),',',num2str(w2(m,n),2),')'), ...
                'horiz','center', ...
                'vert','bottom');
        end
    end
end
end
function plotCS(MD, MC, w1, w2, varargin)

% tamanho da matriz
[M,N] = size(MD);

% valor mínimo e índice
linhaMin = min(MD);
valorMin = min(linhaMin);
%disp(valorMin);
[LMDMin, CMDMin] = find(MD==valorMin); % Linha e Coluna do menor valor na Matriz de Distância

% coordenada mínima de MC
vetorMin = cell2mat(MC(LMDMin,CMDMin)); % vetor com as respectivas coordenada do valor mínimo
%disp(vetorMin);
valorCW1Min = vetorMin(1); % valor da coordenada do menor valor de w1
valorCW2Min = vetorMin(2); % valor da coordenada do menor valor de w2

% valor máximo e índice de MD
linhaMax = max(MD);
valorMax = max(linhaMax);
%disp(valorMax);
[LMDMax, CMDMax] = find(MD==valorMax); % Linha e Coluna do maior valor na Matriz de Distância

% coordenada máxima de MC
vetorMax = cell2mat(MC(LMDMax,CMDMax));
%disp(vetorMax);
valorCW1Max = vetorMax(1); % valor da coordenada do maior valor de w1
valorCW2Max = vetorMax(2); % valor da coordenada do maior valor de w2

% configurações de visualização da matriz
figure(6)
hold on
plot(w1,w2,'r','linewidth',2)
plot(w1',w2','r','linewidth',2)
plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
hold off

% configurações do plot
if strcmp(varargin, 'annotation')
    % x-label, y-label, x-ticks, y-ticks, title
    set(gca, 'FontName', 'Helvetica', 'FontSize', 12)
    xlim([-1 1]);
    ylim([-1 1]);
%     xlabel('w1')
%     ylabel('w2')
    title('Coordenadas Sinápticas')
end

figure(6)
for m = 1:M
    for n = 1:N
        if((w1(m,n)==valorCW1Min) && (w2(m,n)==valorCW2Min))
            text(w1(m,n),w2(m,n),strcat('(',num2str(w1(m,n),2),',',num2str(w2(m,n),2),')'), ...
                'horiz','center', ...
                'vert','bottom', ...
                'FontWeight','bold', ...
                'EdgeColor','black')
        elseif((w1(m,n)==valorCW1Max) && (w2(m,n)==valorCW2Max))
            text(w1(m,n),w2(m,n),strcat('(',num2str(w1(m,n),2),',',num2str(w2(m,n),2),')'), ...
                'horiz','center', ...
                'vert','bottom',...
                'FontWeight','bold', ...
                'EdgeColor','black');
        else
            text(w1(m,n),w2(m,n),strcat('(',num2str(w1(m,n),2),',',num2str(w2(m,n),2),')'), ...
                'horiz','center', ...
                'vert','bottom');
        end
    end
end
end
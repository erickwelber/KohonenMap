function plotMC(MC, MD, varargin)

% dimensão da matriz
[M,N] = size(MC);

% valor mínimo
linhaMin = min(MD);
valorMin = min(linhaMin); % valor mínimo da Matriz de Distâncias
% valor máximo
linhaMax = max(MD);
valorMax = max(linhaMax); % valor máximo da Matriz de Distâncias

% configurações de visualização da matriz
figure(13);
himg = imagesc(MD(:,:));
colormap autumn
% colorbar
grid on

% configurações do plot
if strcmp(varargin, 'annotation')
%     x-label, y-label, x-ticks, y-ticks, title
%     set(gca, 'FontName', 'Helvetica', 'FontSize', 12)
    xlabel('Número da Coluna')
    ylabel('Número da Linha')
    title('Matriz de Coordenadas')
end

for m = 1:M
    for n = 1:N
        if(MD(m,n) == valorMin)
            text(n, m, strcat('(',num2str(MC{m,n},2),')'), ...
                'FontName', 'Helvetica', ...
                'Color','yellow',...
                'FontSize', round(6 + 50./sqrt(M.*N)), ...
                'HorizontalAlignment', 'center', ...
                'Rotation', 0, ...
                'FontWeight','bold', ...
                'EdgeColor','black')
        elseif(MD(m,n) == valorMax)
            text(n, m, strcat('(',num2str(MC{m,n},2),')'), ...
                'FontName', 'Helvetica', ...
                'Color','red',...
                'FontSize', round(6 + 50./sqrt(M.*N)), ...
                'HorizontalAlignment', 'center', ...
                'Rotation', 0, ...
                'FontWeight','bold', ...
                'EdgeColor','black')
        else
            text(n, m, strcat('(',num2str(MC{m,n},2),')'), ...
                'FontName', 'Helvetica', ...
                'FontSize', round(6 + 50./sqrt(M.*N)), ...
                'HorizontalAlignment', 'center', ...
                'Rotation', 0)
        end
    end
end
    % set the datatip UpdateFcn
    cursorMode = datacursormode(gcf);
    set(cursorMode, 'UpdateFcn', {@datatiptxt, himg})
end

function text_to_display = datatiptxt(~, hDatatip, himg)
% determine the current datatip position
pos = get(hDatatip, 'Position');
% form the datatip label
text_to_display = {['Row: ' num2str(pos(2))], ...
    ['Column: ' num2str(pos(1))], ...
    ['Value: ' num2str(himg.CData(pos(2), pos(1)))]};
end
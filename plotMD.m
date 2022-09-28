function plotMD(A, varargin)
% function: matvisual(A, varargin)
% A - m-by-n matrix to be visualized with up to 3 pages
% varargin - type 'annotation' in the place of varargin if one want to
%            annotate the matrix plot (x-label, y-label, etc.)
% check the input
if ~isreal(A) || isempty(A) || ischar(A) || ndims(A) > 3
    errordlg('The data array is not suitable for visualization!', ...
        'Error!', 'modal')
    return
end

% determine the matrix size
[M, N, P] = size(A);

% valor mínimo
linhaMin = min(A);
valorMin = min(linhaMin); % valor mínimo da Matriz de Distâncias
% disp(valorMin);
% valor máximo
linhaMax = max(A);
valorMax = max(linhaMax); % valor máximo da Matriz de Distâncias
% disp(valorMax);

% loop through the matrix pages
for p = 1:P
    
    % prepare page-by-page visualization
    if P > 1, subplot(1, P, p), end
    
    % visualize the matrix
    figure(4)
    himg = imagesc(A(:, :, p));
    colormap autumn
    %     colormap winter;
    %     colormap parula;
    colorbar;
    grid on
    
    % annotation
    if strcmp(varargin, 'annotation')
        % x-label, y-label, x-ticks, y-ticks, title
%         set(gca, 'FontName', 'Helvetica', 'FontSize', 12)
        xlabel('Número da Coluna')
        ylabel('Número da Linha')
        title('Matriz de Distâncias')
        if P > 1, title(['Matrix page ' num2str(p)]), end
        if M <= 50, set(gca, 'YTick', 1:M), end
        if N <= 50, set(gca, 'XTick', 1:N), end
    end
    % values labeling
    for m = 1:M
        for n = 1:N
            if(A(m,n) == valorMin)
                text(n, m, num2str(A(m, n, p), 3), ...
                    'FontName', 'Helvetica', ...
                    'Color','yellow',...
                    'FontSize', round(6 + 50./sqrt(M.*N)), ...
                    'HorizontalAlignment', 'center', ...
                    'Rotation', 0, ...
                    'FontWeight','bold', ...
                    'EdgeColor','black')
            elseif(A(m,n) == valorMax)
                text(n, m, num2str(A(m, n, p), 3), ...
                    'FontName', 'Helvetica', ...
                    'Color','red',...
                    'FontSize', round(6 + 50./sqrt(M.*N)), ...
                    'HorizontalAlignment', 'center', ...
                    'Rotation', 0, ...
                    'FontWeight','bold', ...
                    'EdgeColor','black')
            else
                text(n, m, num2str(A(m, n, p), 3), ...
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
end

function text_to_display = datatiptxt(~, hDatatip, himg)
% determine the current datatip position
pos = get(hDatatip, 'Position');
% form the datatip label
text_to_display = {['Row: ' num2str(pos(2))], ...
    ['Column: ' num2str(pos(1))], ...
    ['Value: ' num2str(himg.CData(pos(2), pos(1)))]};
end

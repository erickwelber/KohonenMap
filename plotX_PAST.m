function plotX_PAST(X_PAST)

% figure(200)
% hold on
% title('Memória de Curto Prazo')
[M,N]=size(X_PAST);
legendInfo={};
for i=1:N
    if(i>=0 && i<=4)
        figure(200)
        hold on
        title('Memória de Curto Prazo')
        plot(300/length(X_PAST(:,i)).*(1:1:length(X_PAST(:,i))),X_PAST(:,i));
        grid on
        legendInfo{end+1} = ['x_j,j= ' num2str(i)];
    end
    legend(legendInfo,'Location','northwest','NumColumns',2)

    if(i>=5 && i<=8)
        figure(201)
        hold on
        title('Memória de Curto Prazo')
        plot(300/length(X_PAST(:,i)).*(1:1:length(X_PAST(:,i))),X_PAST(:,i));
        grid on
        legendInfo{end+1} = ['x_j,j= ' num2str(i)];
    end
    legend(legendInfo,'Location','northwest','NumColumns',2)

    if(i>=9 && i<=12)
        figure(202)
        hold on
        title('Memória de Curto Prazo')
        plot(300/length(X_PAST(:,i)).*(1:1:length(X_PAST(:,i))),X_PAST(:,i));
        grid on
        legendInfo{end+1} = ['x_j,j= ' num2str(i)];
    end
    legend(legendInfo,'Location','northwest','NumColumns',2)

    if(i>=13 && i<=16)
        figure(203)
        hold on
        title('Memória de Curto Prazo')
        plot(300/length(X_PAST(:,i)).*(1:1:length(X_PAST(:,i))),X_PAST(:,i));
        grid on
        legendInfo{end+1} = ['x_j,j= ' num2str(i)];
    end
    legend(legendInfo,'Location','northwest','NumColumns',2)

%     if(i>=201 && i<=250)
%         figure(204)
%         hold on
%         title('Memória de Curto Prazo')
%         plot(300/length(X_PAST(i,:)).*(1:1:length(X_PAST(i,:))),X_PAST);
%         grid on
% %         legendInfo{end+1} = ['Plot ' num2str(i)];
%     end
%     legend(legendInfo,'Location','northwest','NumColumns',2)

%     if(i>=251 && i<=300)
%         figure(205)
%         hold on
%         title('Memória de Curto Prazo')
%         plot(300/length(X_PAST(i,:)).*(1:1:length(X_PAST(i,:))),X_PAST);
%         grid on
% %         legendInfo{end+1} = ['Plot ' num2str(i)];
%     end
% %     legend(legendInfo,'Location','northwest','NumColumns',2)
    
end

hold off

end
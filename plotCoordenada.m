%% Limpa tudo e fecha tudo
clc
clear all
close all

% x = [3,2];
% y = [15,12];
% line(x,y,'LineStyle','--','LineWidth',2);

Mcell = {{[34,42] [82,67], [23,16], [3,2], [15,12]}};
[M,N] = size(Mcell{1,1});
Mcell;

%% Plot
figure(1);
for i=1:N
    if(i+1<=N)
        plot(Mcell{1,1}{1,i}(1),Mcell{1,1}{1,i}(2),'bo','LineWidth',2);
        x = [Mcell{1,1}{1,i}(1),Mcell{1,1}{1,i+1}(1)];
        y = [Mcell{1,1}{1,i}(2),Mcell{1,1}{1,i+1}(2)];
        line(x,y,'LineStyle','--','LineWidth',2);
    else
        plot(Mcell{1,1}{1,i}(1),Mcell{1,1}{1,i}(2),'bo','LineWidth',2);
    end
    hold on
end



% x = (1:5)';
% y = rand(5,1);
% bar(x,y)
% %# show X and Y coordinates
% text(x,y,strcat('(',num2str(x),',',num2str(y,2),')'),'horiz','center','vert','bottom')

%%
clc;
clear all;
close all;

% %%
% w1 = rand(10);
% w2 = rand(10);
% 
% %%
% figure(1)
% plot(w1,w2,'r','linewidth',2)
% hold on
% plot(w1',w2','r','linewidth',2)
% hold on
% for m = 1:10
%     for n = 1:10
%         figure(1)
%         hold on
%         plot(w1(m,n),w2(m,n),'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y');
%         text(w1(m,n),w2(m,n),strcat('(',num2str(w1(m,n),2),',',num2str(w2(m,n),2),')'),'horiz','center','vert','bottom');
%     end
% end

number_of_inputs = 1000;



for n=1:number_of_inputs
    x1(n) = sin(n);
    x2(n) = cos(n);
end

% t = linspace(0,2*pi);

plot(x1,x2,'ob')
% axis([-1 1 -1 1]);        

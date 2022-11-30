clc;
clear all;
close all;

% N=2;
% 
% MATRIZ = cell(N^2-1,N^2-1);
% 
% w1=[0.0537,0.1190;0.0028,0.0176];
% w2=[-0.0665,0.1720;0.0768,-0.0690];
% 
% [m,n] = size(w1);
% linha=0;
% for i=1:m
%     for j=1:n
%         linha = linha + 1;
%         array{linha} = [w1(i,j),w2(i,j)];
%     end
% end
% 
% %subarray = array(2:end);
% 
% %[m,n] = size(MATRIZ);
% 
% [lin,col] = size(array);
% 
% for i=1:col
%     if(i<col)
%         n_indice = array{i}
%         for j=i+1:col
%             n_atual = array{j}
%             if(n_indice ~= n_atual)
%                 if(isempty(MATRIZ{i,j-i}))
%                     MATRIZ{i,j-i} = {[n_indice,n_atual norm(n_indice-n_atual)]};
%                     
%                 elseif([MATRIZ{i,j-i}{end}(1),MATRIZ{i,j-i}{end}(2)] ~= [n_indice] | [MATRIZ{i,j-i}{end}(3),MATRIZ{i,j-i}{end}(4)] ~= [n_atual])
%                     MATRIZ{i,j-i}{end+1} = [[n_indice,n_atual],norm(n_indice-n_atual)];
%                 end
%                 
%             end
%         end
%     end
% end
% 
% w1=[0.0537,0.1190;0.0028,0.1220];
% w2=[-0.0665,0.1720;0.0768,0.1719];
% 
% [m,n] = size(w1);
% linha=0;
% for i=1:m
%     for j=1:n
%         linha = linha + 1;
%         array{linha} = [w1(i,j),w2(i,j)];
%     end
% end
% 
% for i=1:col
%     if(i<col)
%         n_indice = array{i}
%         for j=i+1:col
%             n_atual = array{j}
%             if(n_indice ~= n_atual)
%                 if(isempty(MATRIZ{i,j-i}))
%                     MATRIZ{i,j-i} = {[n_indice,n_atual norm(n_indice-n_atual)]};
%                     
%                 elseif([MATRIZ{i,j-i}{end}(1),MATRIZ{i,j-i}{end}(2)] ~= [n_indice] | [MATRIZ{i,j-i}{end}(3),MATRIZ{i,j-i}{end}(4)] ~= [n_atual])
%                     MATRIZ{i,j-i}{end+1} = [[n_indice,n_atual],norm(n_indice-n_atual)];
%                 end
%                 
%             end
%         end
%     end
% end


% v=[0.0509373363964722,0.0359405353707350,-0.0674776529610739,-0.000327189603571407];
% L=length(v);
% x=sort(interp1(1:L,v,linspace(1,L,100)));

v=[0.0509373363964722,0.0359405353707350,-0.0674776529610739,-0.000327189603571407];
x=-0.02;
[~,index]=min(abs(x-v))

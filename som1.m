%% Script: som1.m
%
% Objetivo:
%
%   Executar a fase de treinamento do Algoritmo de Kohonen
%   Exibir as distâncias percorridas pelos neurônios (Matriz de Distâncias)
%   Exibir as coordenadas sinápticas correspondente à matriz de distância
%   Exibir as trajetórias das coordenadas sinápticas (Menor e Maior)*
%   Calc. a posição (DELTA) de cada neurônio (em cada iteração)*
%       - calc. a velocidade de cada neurônio
%       - visualização das "velocidades" (plot)
%
% Revisão:
%
%   Data            Autor           Descrição
%   ====            =====           =========
%   05/07/2022      Erick           Representação do Mapa em Intervalos de Iterações
%   06/07/2022      Erick           Criação da Matriz de Distâncias (MD)
%   09/07/2022      Erick           Criação da Matriz de Coordenadas (MC)
%   12/07/2022      Erick           Visualização da Matriz de Distâncias
%   12/07/2022      Erick           Visualização da Matriz Coordenadas sinápticas
%   12/07/2022      Erick           Visualização das Coordenadas sinápticas
%
% Referência
% ammar al jodah (2022). self organizing map Kohonen Neural Network (https://www.mathworks.com/matlabcentral/fileexchange/46481-self-organizing-map-kohonen-neural-network), MATLAB Central File Exchange. Retrieved July 13, 2022.
% Hristo Zhivomirov (2022). Matrix Visualization with Matlab (https://www.mathworks.com/matlabcentral/fileexchange/69594-matrix-visualization-with-matlab), MATLAB Central File Exchange. Retrieved July 13, 2022.

%% Apaga tudo, limpa tudo e fecha tudo
clc;
clear all;
close all;

%% Parâmetros iniciais
% parmaters
number_of_inputs=1000; % 1000
N=10; % N^2 points 10

upper_bound_x=1;
lower_bound_x=-1;
upper_bound_m=0.1;
lower_bound_m=-0.1;

alpha_initial=0.1;
sigma_initial=0.05;
neighbour_radius_initial=N/2;

T=300; % number of iterations 300
t=1;

%% Populando os vetores de entrada
% initiate input and neural field
for i=1:number_of_inputs
    x1(i)=rand*(upper_bound_x-lower_bound_x)+lower_bound_x;
    x2(i)=rand*(upper_bound_x-lower_bound_x)+lower_bound_x;
    %     plot(x1,x2,'ob')
    %     hold on
end

%% populando os vetores de pesos
% initiate weight neural
for j1=1:N
    for j2=1:N
        w1(j1,j2)=rand*(upper_bound_m-lower_bound_m)+lower_bound_m;
        w2(j1,j2)=rand*(upper_bound_m-lower_bound_m)+lower_bound_m;
        %         plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        %         plot(w1,w2,'r','linewidth',1)
        %         plot(w1',w2','r','linewidth',1)
        %         hold on
    end
end

%% Matriz de Distâncias
MD = zeros(N,N);

%% Matriz de Coordenadas
MC = cell(10);

%% Visualização inicial [instante t=0]
% initial figures
figure(1)
% input points in the x1,x2 map
plot(x1,x2,'ob')
hold on
% neural field points in the x1,x2 map
plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
plot(w1,w2,'r','linewidth',1)
plot(w1',w2','r','linewidth',1)
hold off
title('t=0');
drawnow

%% start traning
while (t<=T)
    
    % update parameters
    alpha=alpha_initial*(1-t/T);
    sigma=sigma_initial*(1-t/T);
    max_neighbour_radius=round(neighbour_radius_initial*(1-t/T));
    
    % loop over all the input values
    for i=1:number_of_inputs % took one input : a 2D vector
        
        % find minimum distance neural unit (winner)
        e_norm=(x1(i)-w1).^2+(x2(i)-w2).^2; % error distance for each neural node (output error matrix)
        
        minj1=1;minj2=1;
        min_norm=e_norm(minj1,minj2); % select first element in matrix
        
        for j1=1:N
            for j2=1:N
                if e_norm(j1,j2)<min_norm
                    min_norm=e_norm(j1,j2);
                    minj1=j1;
                    minj2=j2;
                end
            end
        end
        
        % winner coordinates
        j1_c= minj1;
        j2_c= minj2;
        
        % coordenada do neurônio vencedor antes da atualização [vencedor]
        cwn1 = [w1(j1_c,j2_c),w2(j1_c,j2_c)]; % coordinate winner neuron 1 [cwn1]
        
        % update the winning neuron
        e_factor = exp(-((j1_c-j1_c).^2+(j2_c-j1_c).^2)/2*sigma);
        w1(j1_c,j2_c)=w1(j1_c,j2_c) + alpha * (x1(i) - w1(j1_c,j2_c));
        w2(j1_c,j2_c)=w2(j1_c,j2_c) + alpha * (x2(i) - w2(j1_c,j2_c));
        
        % coordenada do neurônio vencedor depois da atualização [vencedor]
        cwn2 = [w1(j1_c,j2_c),w2(j1_c,j2_c)]; % coordinate winner neuron 2 [cwn2]
        
        % atualiza a coordenada do neurônio [vencedor]
        MC(j1_c,j2_c) = {(cwn2)}; % Matriz de Coordenadas [MC]
        
        % incrementa e atualiza a distância entre os neurônios [vencedor]
        MD(j1_c,j2_c) = MD(j1_c,j2_c) + norm(cwn1-cwn2); % Matriz de Distâncias [MD]
        
        % update the neighbour neurons
        for neighbour_radius=1:1:max_neighbour_radius
            jj1=j1_c - neighbour_radius;
            jj2=j2_c;
            if (jj1>=1) % to stay in the matrix
                % coordenada do neurônio vizinho antes da atualização
                cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
                
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % atualiza a coordenada do neurônio [vizinha]
                MC(jj1,jj2) = {(cnn2)}; % Matriz de Coordenadas [MC]
                
                % incrementa e atualiza a distância entre os neurônios [vizinho]
                MD(jj1,jj2) = MD(jj1,jj2) + norm(cnn1-cnn2); % Matriz de Distâncias [MD]
            end
            jj1=j1_c + neighbour_radius;
            jj2=j2_c;
            if (jj1<=N) % to stay in the matrix
                % coordenada do neurônio vizinho antes da atualização
                cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
                
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % atualiza a coordenada do neurônio [vizinha]
                MC(jj1,jj2) = {(cnn2)}; % Matriz de Coordenadas [MC]
                
                % incrementa e atualiza a distância entre os neurônios [vizinho]
                MD(jj1,jj2) = MD(jj1,jj2) + norm(cnn1-cnn2); % Matriz de Distâncias [MD]
            end
            jj1=j1_c;
            jj2=j2_c - neighbour_radius;
            if (jj2>=1) % to stay in the matrix
                % coordenada do neurônio vizinho antes da atualização
                cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
                
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % atualiza a coordenada do neurônio [vizinha]
                MC(jj1,jj2) = {(cnn2)}; % Matriz de Coordenadas [MC]
                
                % incrementa e atualiza a distância entre os neurônios [vizinho]
                MD(jj1,jj2) = MD(jj1,jj2) + norm(cnn1-cnn2); % Matriz de Distâncias [MD]
            end
            jj1=j1_c;
            jj2=j2_c + neighbour_radius;
            if (jj2<=N) % to stay in the matrix
                % coordenada do neurônio vizinho antes da atualização
                cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
                
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % atualiza a coordenada do neurônio [vizinha]
                MC(jj1,jj2) = {(cnn2)}; % Matriz de Coordenadas [MC]
                
                % incrementa e atualiza a distância entre os neurônios [vizinho]
                MD(jj1,jj2) = MD(jj1,jj2) + norm(cnn1-cnn2); % Matriz de Distâncias [MD]
            end
        end
        
        %         figure(1)
        %         % input points in the x1,x2 map
        %         plot(x1,x2,'ob')
        %         hold on
        %         % neural field points in the x1,x2 map
        %         plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        %         plot(w1,w2,'r','linewidth',1)
        %         plot(w1',w2','r','linewidth',1)
        %         hold off
        %         title(num2str(i));
        %         drawnow
    end
    
    %% plot em intervalo de iterações
    if(t==5)
        figure(3)
        subplot(3,3,1)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==10)
        figure(3)
        subplot(3,3,2)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==20)
        figure(3)
        subplot(3,3,3)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==50)
        figure(3)
        subplot(3,3,4)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==100)
        figure(3)
        subplot(3,3,5)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==150)
        figure(3)
        subplot(3,3,6)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==200)
        figure(3)
        subplot(3,3,7)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==250)
        figure(3)
        subplot(3,3,8)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end
    if(t==300)
        figure(3)
        subplot(3,3,9)
        plot(x1,x2,'ob')
        hold on
        plot(w1,w2,'r','linewidth',2)
        plot(w1',w2','r','linewidth',2)
        plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
        hold off
        title(['t=' num2str(t)]);
        drawnow
    end

    figure(2)
    plot(x1,x2,'ob')
    hold on
    plot(w1,w2,'r','linewidth',2)
    plot(w1',w2','r','linewidth',2)
    plot(w1,w2,'ro','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
    hold off
    title(['t=' num2str(t)]);
    drawnow
    t=t+1;
end

%% visualização da Matriz de Distâncias no plot
plotMD(MD,'annotation')

%% visualização da Matriz de Coordenadas no plot
plotMC(MC,MD,'annotation')

%% visualização das Coordenadas Sinápticas no plot
plotCS(MD, MC, w1, w2, 'annotation')
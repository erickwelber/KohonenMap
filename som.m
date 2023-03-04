%% Script: som1.m
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
number_of_inputs=20; % 1000
N=4; % N^2 points 10

upper_bound_x=1;
lower_bound_x=-1;
upper_bound_m=0.1;
lower_bound_m=-0.1;

alpha_initial=0.1;
sigma_initial=0.05;
neighbour_radius_initial=N/2;

T=300; % número de interações 300
t=1;

% Configuração da memória de curto prazo
xj=0;
xj_past=0;
yi=zeros(1,N^2);
Bj=rand(1,N^2);
aj=0.5;
Dij=rand(N^2,N^2);
x_past=zeros(1,N^2);
MD_temp=zeros(1,N^2);

% parâmetros do plot Gaussiano
max_radius=round(neighbour_radius_initial*(1-t/T));
Amplitude = zeros(N,N);

%% Populando os vetores de entrada CÍRCULO VAZADO
% for n=1:number_of_inputs
%     x1(n) = cos(n);
%     x2(n) = sin(n);
% end

%% Populando os vetores de entrada TRIÂNGULO PREENCHIDO
% x = [-1 1]; % //triangle base
% y = [-1 1]; % //triangle height
% 
% points = rand(number_of_inputs,2); %// sample uniformly in unit square
% ind = points(:,2)>points(:,1); %// points to be unfolded 
% points(ind,:) = [2-points(ind,2) points(ind,1)]; %// unfold them
% points(:,1) = x(1) + (x(2)-x(1))/2*points(:,1); %// stretch x as needed
% points(:,2) = y(1) + (y(2)-y(1))*points(:,2); %// stretch y as needed
% x1 = points(:,1);
% x2 = points(:,2);

%% Populando os vetores de entrada CÍRCULO PREENCHIDO
% for n=1:number_of_inputs
%     
%     r(n) = sqrt(rand(1,1));
%     % and theta as before:
%     theta(n)=2*pi*rand(1,1);
%    
%     % convert to cartesian
%     x1(n)=r(n)*cos(theta(n));
%     x2(n)=r(n)*sin(theta(n));
%     
% end

%% Populando os vetores de entrada RETÂNGULO PREENCHIDO
%initiate input and neural field
for i=1:number_of_inputs
    x1(i)=rand*(upper_bound_x-lower_bound_x)+lower_bound_x;
    x2(i)=rand*(upper_bound_x-lower_bound_x)+lower_bound_x;
end

%% Populando os vetores de pesos
% initiate weight neural
for j1=1:N
    for j2=1:N
        w1(j1,j2)=rand*(upper_bound_m-lower_bound_m)+lower_bound_m;
        w2(j1,j2)=rand*(upper_bound_m-lower_bound_m)+lower_bound_m;
    end
end

%% Matriz de Distâncias
MD = zeros(N,N);

%% Matriz de Coordenadas
MC = cell(N);

%% Matriz de Trajetórias
MT = cell(N);

%% Matriz de Distâncias entre Neurônios em Cada Atualização
MDEN = cell(N^2-1,N^2-1);

%% Matriz de Distâncias entre TODOS os Neurônios em Cada Atualização
MDETN = cell(N,N);

%% Visualização inicial [instante t=0]
% initial figures
figure(1)
hold on
title('t=0');
% input points in the x1,x2 map
plot(x1,x2,'ob')
% neural field points in the x1,x2 map
plot(w1,w2,'r','linewidth',1)
plot(w1',w2','r','linewidth',1)
plot(w1,w2,'yo','MarkerSize',6,'MarkerEdgeColor','r','MarkerFaceColor','y')
hold off

%% Duração inicial
t_inicial = cputime;

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
        
        % Acumulação das coordenadas vencedoras antes da atualização
        if(isempty(MT{j1_c,j2_c}))
            MT{j1_c,j2_c} = {cwn1};
        elseif(MT{j1_c,j2_c}{end} ~= cwn1)
            MT{j1_c,j2_c}{end+1} = (cwn1);
        end
        
        % parameters signal Gaussian winner
        [signal_winner, value_neighbour] = gaussianW(w1,w2,j1_c,j2_c,max_neighbour_radius,N);
        Amplitude(j1_c,j2_c) = signal_winner;

        % update the winning neuron
        e_factor = exp(-((j1_c-j1_c).^2+(j2_c-j1_c).^2)/2*sigma);
        w1(j1_c,j2_c)=w1(j1_c,j2_c) + alpha * (x1(i) - w1(j1_c,j2_c)) + signal_winner;
        w2(j1_c,j2_c)=w2(j1_c,j2_c) + alpha * (x2(i) - w2(j1_c,j2_c)) + signal_winner;
        
        % coordenada do neurônio vencedor depois da atualização [vencedor]
        cwn2 = [w1(j1_c,j2_c),w2(j1_c,j2_c)]; % coordinate winner neuron 2 [cwn2]
        
        % Acumulação das coordenadas vencedoras depois da atualização
        if(isempty(MT{j1_c,j2_c}))
            MT{j1_c,j2_c} = {cwn2};
        elseif(MT{j1_c,j2_c}{end} ~= cwn2)
            MT{j1_c,j2_c}{end+1} = (cwn2);            
        end
        
        % Matriz de Coordenadas entre neurônios depois da atualização
        [m,n] = size(w1);
        linha = 0;
        for i=1:m
            for j=1:n
                linha = linha + 1;
                array{linha} = [w1(i,j),w2(i,j)];
            end
        end
        [~,col] = size(array);
        for i=1:col
            if(i<col)
                n_indice = array{i};
                for j=i+1:col
                    n_atual = array{j};
                    if(n_indice ~= n_atual)
                        if(isempty(MDEN{i,j-i}))
                            MDEN{i,j-i} = {[n_indice,n_atual norm(n_indice-n_atual)]};
                        elseif([MDEN{i,j-i}{end}(1),MDEN{i,j-i}{end}(2)] ~= [n_indice] | [MDEN{i,j-i}{end}(3),MDEN{i,j-i}{end}(4)] ~= [n_atual])
                            MDEN{i,j-i}{end+1} = [[n_indice,n_atual],norm(n_indice-n_atual)];
                        end
                    end
                end
            end
        end

        %------------------------------------------------------------------
        % matriz de coordenada entre todos os neurônios
        [m,n] = size(w1);
        linha = 0;
        for i=1:m
            for j=1:n
                linha = linha + 1;
                MP{linha} = [w1(i,j),w2(i,j)];
            end
        end
        [~,col] = size(MP);
        for i=1:col
            n_indice = MP{i};
            for j=1:col
                n_atual = MP{j};
                MDETN{i,j} = {[n_indice,n_atual,norm(n_indice-n_atual)]};
            end
        end
        % obtendo índice de yi
        TAM = length(MP);
        for ii=1:TAM
            if(cwn2==MP{ii})
                j_index=ii;
                break
            end
        end
        yi(j_index)= Amplitude(j1_c,j2_c);
        % obtendo yi
%         for ii=1:TAM
%             yi(ii) = MDETN{indice_yi,ii}{1}(5);
%         end
        %------------------------------------------------------------------

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
                
                % Acumulação das coordenadas vizinhas antes da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn1};
                elseif(MT{jj1,jj2}{end} ~= cnn1)
                    MT{jj1,jj2}{end+1} = (cnn1);
                end
                
                % parameters signal Gaussian neighbour
                [signal_neighbour] = gaussianN(w1,w2,jj1,jj2);
                Amplitude(jj1,jj2) = value_neighbour(jj1,jj2) + signal_neighbour;

                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % Acumulação das coordenadas vizinhas depois da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn2};
                elseif(MT{jj1,jj2}{end} ~= cnn2)
                    MT{jj1,jj2}{end+1} = (cnn2);
                end
                
                % Matriz de Coordenadas entre neurônios depois da atualização
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        array{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [lin,col] = size(array);
                for i=1:col
                    if(i<col)
                        n_indice = array{i};
                        for j=i+1:col
                            n_atual = array{j};
                            if(n_indice ~= n_atual)
                                if(isempty(MDEN{i,j-i}))
                                    MDEN{i,j-i} = {[n_indice,n_atual norm(n_indice-n_atual)]};
                                elseif([MDEN{i,j-i}{end}(1),MDEN{i,j-i}{end}(2)] ~= [n_indice] | [MDEN{i,j-i}{end}(3),MDEN{i,j-i}{end}(4)] ~= [n_atual])
                                    MDEN{i,j-i}{end+1} = [[n_indice,n_atual],norm(n_indice-n_atual)];
                                end
                            end
                        end
                    end
                end

                %----------------------------------------------------------
                % matriz de coordenada entre todos os neurônios
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        MP{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [~,col] = size(MP);
                for i=1:col
                    n_indice = MP{i};
                    for j=1:col
                        n_atual = MP{j};
                        MDETN{i,j} = {[n_indice,n_atual,norm(n_indice-n_atual)]};
                    end
                end
                % obtendo índice de yi
                TAM = length(MP);
                for ii=1:TAM
                    if(cnn2==MP{ii})
                        indice_yi=ii;
                        break
                    end
                end
                yi(indice_yi)= Amplitude(jj1,jj2);
                MD_temp(indice_yi) = MDETN{j_index,indice_yi}{1}(5);
                % obtendo yi
%                 for ii=1:TAM
%                     yi(ii) = MDETN{indice_yi,ii}{1}(5);
%                 end
                %----------------------------------------------------------
                
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
                
                % Acumulação das coordenadas vizinhas antes da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn1};
                elseif(MT{jj1,jj2}{end} ~= cnn1)
                    MT{jj1,jj2}{end+1} = (cnn1);
                end
                
                % parameters signal Gaussian neighbour
                [signal_neighbour] = gaussianN(w1,w2,jj1,jj2);
                Amplitude(jj1,jj2) = value_neighbour(jj1,jj2) + signal_neighbour;
                
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % Acumulação das coordenadas vizinhas depois da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn2};
                elseif(MT{jj1,jj2}{end} ~= cnn2)
                    MT{jj1,jj2}{end+1} = (cnn2);
                end
                
                % Matriz de Coordenadas entre neurônios depois da atualização
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        array{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [lin,col] = size(array);
                for i=1:col
                    if(i<col)
                        n_indice = array{i};
                        for j=i+1:col
                            n_atual = array{j};
                            if(n_indice ~= n_atual)
                                if(isempty(MDEN{i,j-i}))
                                    MDEN{i,j-i} = {[n_indice,n_atual norm(n_indice-n_atual)]};
                                elseif([MDEN{i,j-i}{end}(1),MDEN{i,j-i}{end}(2)] ~= [n_indice] | [MDEN{i,j-i}{end}(3),MDEN{i,j-i}{end}(4)] ~= [n_atual])
                                    MDEN{i,j-i}{end+1} = [[n_indice,n_atual],norm(n_indice-n_atual)];
                                end
                            end
                        end
                    end
                end

                %----------------------------------------------------------
                % matriz de coordenada entre todos os neurônios
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        MP{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [~,col] = size(MP);
                for i=1:col
                    n_indice = MP{i};
                    for j=1:col
                        n_atual = MP{j};
                        MDETN{i,j} = {[n_indice,n_atual,norm(n_indice-n_atual)]};
                    end
                end
                % obtendo índice de yi
                TAM = length(MP);
                for ii=1:TAM
                    if(cnn2==MP{ii})
                        indice_yi=ii;
                        break
                    end
                end
                yi(indice_yi)= Amplitude(jj1,jj2);
                MD_temp(indice_yi) = MDETN{j_index,indice_yi}{1}(5);
                % obtendo yi
%                 for ii=1:TAM
%                     yi(ii) = MDETN{indice_yi,ii}{1}(5);
%                 end
                %----------------------------------------------------------

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
                
                % Acumulação das coordenadas vizinhas antes da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn1};
                elseif(MT{jj1,jj2}{end} ~= cnn1)
                    MT{jj1,jj2}{end+1} = (cnn1);
                end
                
                % parameters signal Gaussian neighbour
                [signal_neighbour] = gaussianN(w1,w2,jj1,jj2);
                Amplitude(jj1,jj2) = value_neighbour(jj1,jj2) + signal_neighbour;
                
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % Acumulação das coordenadas vizinhas depois da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn2};
                elseif(MT{jj1,jj2}{end} ~= cnn2)
                    MT{jj1,jj2}{end+1} = (cnn2);
                end
                
                % Matriz de Coordenadas entre neurônios depois da atualização
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        array{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [lin,col] = size(array);
                for i=1:col
                    if(i<col)
                        n_indice = array{i};
                        for j=i+1:col
                            n_atual = array{j};
                            if(n_indice ~= n_atual)
                                if(isempty(MDEN{i,j-i}))
                                    MDEN{i,j-i} = {[n_indice,n_atual norm(n_indice-n_atual)]};
                                elseif([MDEN{i,j-i}{end}(1),MDEN{i,j-i}{end}(2)] ~= [n_indice] | [MDEN{i,j-i}{end}(3),MDEN{i,j-i}{end}(4)] ~= [n_atual])
                                    MDEN{i,j-i}{end+1} = [[n_indice,n_atual],norm(n_indice-n_atual)];
                                end
                            end
                        end
                    end
                end

                %----------------------------------------------------------
                % matriz de coordenada entre todos os neurônios
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        MP{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [~,col] = size(MP);
                for i=1:col
                    n_indice = MP{i};
                    for j=1:col
                        n_atual = MP{j};
                        MDETN{i,j} = {[n_indice,n_atual,norm(n_indice-n_atual)]};
                    end
                end
                % obtendo índice de yi
                TAM = length(MP);
                for ii=1:TAM
                    if(cnn2==MP{ii})
                        indice_yi=ii;
                        break
                    end
                end
                yi(indice_yi)= Amplitude(jj1,jj2);
                MD_temp(indice_yi) = MDETN{j_index,indice_yi}{1}(5);
                % obtendo yi
%                 for ii=1:TAM
%                     yi(ii) = MDETN{indice_yi,ii}{1}(5);
%                 end
                %----------------------------------------------------------
                
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
                
                % Acumulação das coordenadas vizinhas antes da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn1};
                elseif(MT{jj1,jj2}{end} ~= cnn1)
                    MT{jj1,jj2}{end+1} = (cnn1);
                end
                
                % parameters signal Gaussian neighbour
                [signal_neighbour] = gaussianN(w1,w2,jj1,jj2);
                Amplitude(jj1,jj2) = value_neighbour(jj1,jj2) + signal_neighbour;
                
                e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
                w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2)) + value_neighbour(jj1,jj2) + signal_neighbour;
                
                % coordenada do neurônio vizinho depois da atualização
                cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
                
                % Acumulação das coordenadas vizinhas depois da atualização
                if(isempty(MT{jj1,jj2}))
                    MT{jj1,jj2} = {cnn2};
                elseif(MT{jj1,jj2}{end} ~= cnn2)
                    MT{jj1,jj2}{end+1} = (cnn2);
                end
                
                % Matriz de Coordenadas entre neurônios depois da atualização
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        array{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [lin,col] = size(array);
                for i=1:col
                    if(i<col)
                        n_indice = array{i};
                        for j=i+1:col
                            n_atual = array{j};
                            if(n_indice ~= n_atual)
                                if(isempty(MDEN{i,j-i}))
                                    MDEN{i,j-i} = {[n_indice,n_atual norm(n_indice-n_atual)]};
                                elseif([MDEN{i,j-i}{end}(1),MDEN{i,j-i}{end}(2)] ~= [n_indice] | [MDEN{i,j-i}{end}(3),MDEN{i,j-i}{end}(4)] ~= [n_atual])
                                    MDEN{i,j-i}{end+1} = [[n_indice,n_atual],norm(n_indice-n_atual)];
                                end
                            end
                        end
                    end
                end

                %----------------------------------------------------------
                % matriz de coordenada entre todos os neurônios
                [m,n] = size(w1);
                linha = 0;
                for i=1:m
                    for j=1:n
                        linha = linha + 1;
                        MP{linha} = [w1(i,j),w2(i,j)];
                    end
                end
                [~,col] = size(MP);
                for i=1:col
                    n_indice = MP{i};
                    for j=1:col
                        n_atual = MP{j};
                        MDETN{i,j} = {[n_indice,n_atual,norm(n_indice-n_atual)]};
                    end
                end
                % obtendo índice de yi
                TAM = length(MP);
                for ii=1:TAM
                    if(cnn2==MP{ii})
                        indice_yi=ii;
                        break
                    end
                end
                yi(indice_yi)= Amplitude(jj1,jj2);
                MD_temp(indice_yi) = MDETN{j_index,indice_yi}{1}(5);
                % obtendo yi
%                 for ii=1:TAM
%                     yi(ii) = MDETN{indice_yi,ii}{1}(5);
%                 end
                %----------------------------------------------------------
                
                % atualiza a coordenada do neurônio [vizinha]
                MC(jj1,jj2) = {(cnn2)}; % Matriz de Coordenadas [MC]
                
                % incrementa e atualiza a distância entre os neurônios [vizinho]
                MD(jj1,jj2) = MD(jj1,jj2) + norm(cnn1-cnn2); % Matriz de Distâncias [MD]
            end
        end

        options = odeset('RelTol',1e-4,'AbsTol',[1e-4]);
        [tint,xj] = ode45(@(t,xj) x_stm(t,xj,MD_temp,yi,Bj(j_index),Dij(:,j_index),aj,x_past,N),[0 1], 0,options);
        yi = zeros(1,N^2);
        x_past(j_index) = xj(length(xj));
    end

% plot em intervalo de iterações
% plotITERA(x1,x2,w1,w2,t);
% plotITERA(x1,x2,w1,w2,t,max_radius,j1_c,j2_c,N,Amplitude);

X_PAST(t,:) = x_past(1,:);

t=t+1;
end

%%Duração final
t_final = cputime;
total = t_final - t_inicial

%% visualização do mapa completo
plotMAP(x1,x2,w1,w2,t);

%% visualização da Matriz de Distâncias no plot
% plotMD(MD, 'annotation');

%% visualização da Matriz de Coordenadas no plot
% plotMC(MC, MD, 'annotation');

%% visualização das Coordenadas Sinápticas no plot
% plotCS(MD, MC, w1, w2);

%% visualização das trajetórias mínimas e máximas
% plotTJ(MT, MD);

%% visualização dos neurônios 1 e 2
% plotMDEN(MDEN);

%% visualização da gaussiana em 3D
% plotGAUSS(x1,x2,w1,w2,max_radius,j1_c,j2_c,N,Amplitude);

% plotar x_past para cada iteração 't'
plotX_PAST(X_PAST);

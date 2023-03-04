function [signal_winner, value_neighbour] = gaussianW(w1,w2,j1_c,j2_c,max_neighbour_radius,N)

% definições preliminares
x = linspace(-1,1,100);
y = linspace(-1,1,100);
[X,Y] = meshgrid(x,y);

% x=reshape(w1,[1,numel(w1)]);
% Lx=length(x);
% xx=sort(interp1(1:Lx,x,linspace(1,Lx,100)));
% y=reshape(w2,[1,numel(w2)]);
% Ly=length(y);
% yy=sort(interp1(1:Ly,y,linspace(1,Ly,100)));
% [X,Y] = meshgrid(xx,yy);

% parâmetros
ux = w1(j1_c,j2_c);
uy = w2(j1_c,j2_c);
sig = 0.5;
Amp = 10*(1/(sqrt(2*pi*sig^2)));

% Gaussian
fx = Amp * exp(-(X-ux).^2/(2*sig^2));
fy = Amp * exp(-(Y-uy).^2/(2*sig^2));

gauss2d = fx.*fy;

cx=X(1,:);
cy=Y(:,1)';

% signal_winner = max(max(gauss2d(:)))
% [~,indexcX] = min(abs(ux-cx));     
% [~,indexcY] = min(abs(uy-cy));
% signal_winner_2 = gauss2d(indexcY,indexcX)

value_neighbour = NaN(N);

% -----
for neighbour_radius=1:1:max_neighbour_radius
    jj1=j1_c - neighbour_radius;
    jj2=j2_c;
    if (jj1>=1) % to stay in the matrix
        nx=w1(jj1,jj2);
        ny=w2(jj1,jj2);
        [~,indexX] = min(abs(nx-cx));      
        [~,indexY] = min(abs(ny-cy));
        value_neighbour(jj1,jj2) = gauss2d(indexY,indexX);
%         value_neighbour = value_neighbour + gauss2d(indexY,indexX);
%         % coordenada do neurônio vizinho antes da atualização
%         cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
%         
%         e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
%         w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
%         w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
%         
%         % coordenada do neurônio vizinho depois da atualização
%         cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
    end
    jj1=j1_c + neighbour_radius;
    jj2=j2_c;
    if (jj1<=N) % to stay in the matrix
        nx=w1(jj1,jj2);
        ny=w2(jj1,jj2);
        [~,indexX] = min(abs(nx-cx));      
        [~,indexY] = min(abs(ny-cy));
        value_neighbour(jj1,jj2) = gauss2d(indexY,indexX);
%         value_neighbour = value_neighbour + gauss2d(indexY,indexX);
%         % coordenada do neurônio vizinho antes da atualização
%         cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
%         
%         e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
%         w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
%         w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
%         
%         % coordenada do neurônio vizinho depois da atualização
%         cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
    end
    jj1=j1_c;
    jj2=j2_c - neighbour_radius;
    if (jj2>=1) % to stay in the matrix
        nx=w1(jj1,jj2);
        ny=w2(jj1,jj2);
        [~,indexX] = min(abs(nx-cx));      
        [~,indexY] = min(abs(ny-cy));
        value_neighbour(jj1,jj2) = gauss2d(indexY,indexX);
%         value_neighbour = value_neighbour + gauss2d(indexY,indexX);
%         % coordenada do neurônio vizinho antes da atualização
%         cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
%         
%         e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
%         w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
%         w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
%         
%         % coordenada do neurônio vizinho depois da atualização
%         cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
    end
    jj1=j1_c;
    jj2=j2_c + neighbour_radius;
    if (jj2<=N) % to stay in the matrix
        nx=w1(jj1,jj2);
        ny=w2(jj1,jj2);
        [~,indexX] = min(abs(nx-cx));      
        [~,indexY] = min(abs(ny-cy));
        value_neighbour(jj1,jj2) = gauss2d(indexY,indexX);
%         value_neighbour = value_neighbour + gauss2d(indexY,indexX);
%         % coordenada do neurônio vizinho antes da atualização
%         cnn1 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 1 [cnn1]
%         
%         e_factor = exp(-((j1_c-jj1).^2+(j2_c-jj2).^2)/2*sigma);
%         w1(jj1,jj2)=w1(jj1,jj2) + alpha * e_factor * (x1(i)-w1(jj1,jj2));
%         w2(jj1,jj2)=w2(jj1,jj2) + alpha * e_factor * (x2(i)-w2(jj1,jj2));
%         
%         % coordenada do neurônio vizinho depois da atualização
%         cnn2 = [w1(jj1,jj2),w2(jj1,jj2)]; % coordenate neighbour neuron 2 [cnn2]
    end
end
% -----

signal_winner = max(max(gauss2d(:)));

% figure(20), clf
% surf(x,y,gauss2d);
% colormap autumn;
% axis square;
% rotate3d on;
% shading interp

% plotMD(gauss2d,'annotation');
end
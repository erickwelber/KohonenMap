function [signal_winner_neighbour] = gaussianN(w1,w2,jj1,jj2)

% definições preliminares
x = linspace(-1,1);
y = linspace(-1,1);
[X,Y] = meshgrid(x,y);

% x=reshape(w1,[1,numel(w1)]);
% Lx=length(x);
% xx=sort(interp1(1:Lx,x,linspace(1,Lx,100)));
% y=reshape(w2,[1,numel(w2)]);
% Ly=length(y);
% yy=sort(interp1(1:Ly,y,linspace(1,Ly,100)));
% [X,Y] = meshgrid(xx,yy);

% parâmetros
ux = w1(jj1,jj2);
uy = w2(jj1,jj2);
sig = 0.5;
Amp = 0.5*(1/(sqrt(2*pi*sig^2)));

% Gaussian
fx = Amp * exp(-(X-ux).^2/(2*sig^2));
fy = Amp * exp(-(Y-uy).^2/(2*sig^2));

gauss2d = fx.*fy;

signal_winner_neighbour = max(max(gauss2d(:)));

% figure(30), clf
% surf(x,y,gauss2d);
% colormap autumn;
% axis square;
% rotate3d on;
% shading interp

% plotMD(gauss2d,'annotation');
end
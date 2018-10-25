
% Implementación de Kalman para seguir la trayectoria de un vehículo

close all;
clc

%% Mido posición con ruido blanco de sigma = 0.01

sigma = 0.01;
C = [eye(3) zeros(3) zeros(3)];
x = [Pos'; Vel'; Acel'];
eta = sigma.*randn(3,length(Pos));

y_pos = C*x + eta; % obtengo medición de la posición con ruido 

x_0_0 = x(1) + sigma.*randn(9,1);
P_0_0 = [eye(3)*sigma.*randn(1) zeros(3) zeros(3);...
    zeros(3) eye(3)*sigma/10.*randn(1) zeros(3);...
    zeros(3) zeros(3) eye(3)*sigma/100.*randn(1)];

h = abs(tiempos(2) - tiempos(1));
A = [eye(3) eye(3)*h eye(3)*h^2;...
    zeros(3) eye(3) eye(3)*h;...
    zeros(3) zeros(3) eye(3)];





%% Trayectoria de un cohete

clear all;
close all;

trayectoria = dlmread('trayectoria.csv',',',1,0);
tiempo = trayectoria(:,1)';
Pos = trayectoria(:,2:4)';
Vel = trayectoria(:,5:7)';
Acel = trayectoria(:,8:10)';
clear trayectoria;

    %% Determinación del sistema

I = eye(3);
O = zeros(3);
h = abs(tiempo(1)-tiempo(2));

x = [Pos; Vel; Acel];

A = [I  I*h  I*(h^2)/2;
     O  I    I*h;
     O  O    I];
 
% Ruido del proceso: PREGUNTAR CRITERIO
B = [I  O  O;
     O  I  O;
     O  O  I];
 
sigma_p = 100; % Ruido blanco 
Q = [I*sigma_p^2    O               O;
     O              I*sigma_p^2     O;
     O              O               I*sigma_p^2];


%% 1) Se mide la posición afectada por ruido blanco (gaussiano) de 10m de desvı́ío estándar

C = [I O O]; % mido posición

sigma_m = 10; % desvío estándar de medición
R = I*sigma_m^2; % Matriz de correlación (diagonal por ruido blanco)

% Genero una matriz que tenga en cada columna una realización de un proceso
% de ruido blanco gaussiano con matriz de covarianza R y media cero.
eta = ruido_blanco_gaussiano(R,length(x)); 

% Genero las muestras con ruido blanco gaussiano R.
y = C*x + eta;

xhat = zeros(size(x));
P = [I  O  O;
     O  I  O;
     O  O  I];
  
for n=1:length(x)-1
    xhat(:,n+1) = A*xhat(:,n); % Predicción
    P = A*P*A' + B*Q*B';       
    K = P*C'*(inv(R+C*P*C'));
    xhat(:,n+1) = xhat(:,n+1) + K*(y(:,n+1) -C*xhat(:,n+1)); % Actualización
    P = (eye(length(A)) - K*C)*P;
end

figure(1); hold all;
title('Medición de la posición afectada por ruido blanco gaussiano de 10m de desvı́ío estándar')
subplot(3,1,1)
plot(tiempo,x(1,:),tiempo,xhat(1,:),tiempo,y(1,:),'go')
xlabel('Tiempo')
ylabel('Px')
grid on;
legend('Proceso real','Proceso estimado','medición')
subplot(3,1,2)
plot(tiempo,x(2,:),tiempo,xhat(2,:),tiempo,y(1,:),'go')
xlabel('Tiempo')
ylabel('Py')
grid on;
legend('Proceso real','Proceso estimado','medición')
subplot(3,1,3)
plot(tiempo,x(3,:),tiempo,xhat(3,:),tiempo,y(1,:),'go')
xlabel('Tiempo')
ylabel('Pz')
grid on;
legend('Proceso real','Proceso estimado','medición')

%% 2) Se mide la posición afectada por ruido blanco (uniforme) de 10m de desvı́ío estándar.






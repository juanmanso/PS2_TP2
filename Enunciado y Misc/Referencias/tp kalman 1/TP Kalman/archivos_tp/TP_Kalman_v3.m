%% Kalman sin sesgo

clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
h = 1; % Período de muestreo.
Abx = Acel(1:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(1:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(1:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz
N = length(Abx); % Cantidad de muestras

% La dimensión del vector $x$ a estimar es n = dim(Pos)+dim(Vel)+dim(Cb), 
% la del ruido $\xi$ del proceso es q = dim(Pos)+dim(Vel)+dim(Cb), y la 
% del ruido $\eta$ de medición es p = dim(Pos)+dim(Vel).
n = size(Preal,2)-1 + size(Vreal,2)-1 + 4;
q = size(Preal,2)-1 + size(Vreal,2)-1 + 4;
p = size(Preal,2)-1 + size(Vreal,2)-1;

% Matrices útiles:
I = eye(size(Preal,2)-1);
O = zeros(size(Preal,2)-1);

% Dinámica del proceso:
A = zeros(N,n,n);
for i=1:N
    A(i,:,:) = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
                0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
                0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
                0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
                0 0 0 0 h             w(i)*h         0             0;
                0 0 0 0 -w(i)*h       h              0             0;
                0 0 0 0 0             0              h             w(i)*h;
                0 0 0 0 0             0              -w(i)*h       h];
end

% Ruido del proceso:
B = [I*h^3/6 O       O       O;
     O       I*h^2/2 O       O;
     O       O       I*h^2/2 O;
     O       O       O       I*h^2/2];
 
s_proceso = .00001; % desvío del ruido del proceso.
Q = eye(q)*s_proceso;

% Medición. Se mide posición y velocidad con el radar y se guarda cada
% medición en una columna del vector y:
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
C = [I O O O;
     O I O O]; 

% Ruido de medición:
sp = 10;
sv = .1;

R = [sp^2*I O;   
     O      sv^2*I];

% Condiciones iniciales:
deltaP = 100;
deltaV = .2;
deltaTheta = 2/9*pi;
deltaC = max(sin(deltaTheta),cos(deltaTheta));
P = [deltaP*I O        O        O;
     O        deltaV*I O        O;
     O        O        deltaC*I O;
     O        O        O        deltaC*I];
x0 = zeros(n,1);

xhat = kalman1(A,B,C,0,0,0,Q,R,x0,P,y); % Estimación por ecuaciones de Riccati
%xhat = kalman2(A,B,C,Q,R,x0,P,y); % Estimación por factorización QR

% Gráficos:
figure(1); %gráfico de la posición
subplot(2,2,1);
plot(Preal(1:100:end-1,1),Preal(1:100:end-1,2),Preal(1:100:end-1,1),xhat(1,:));
title('Posición en x');
grid on;
legend('Real','Estimada');
subplot(2,2,2);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,2)-xhat(1,:)'));
ylim([0 500]);
title('Error de estimación en x');
grid on;
subplot(2,2,3);
plot(Preal(1:100:end-1,1),Preal(1:100:end-1,3),Preal(1:100:end-1,1),xhat(2,:));
title('Posición en y');
grid on;
legend('Real','Estimada');
subplot(2,2,4);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,3)-xhat(2,:)'));
ylim([0 500]);
title('Error de estimación en y');
grid on;


figure(2); %gráfico de la velocidad en x y en y en función del tiempo
subplot(2,2,1);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),xhat(3,:));
title('Velocidad en x');
grid on;
legend('Real','Estimada');
subplot(2,2,2);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,2)-xhat(3,:)'));
ylim([0 5]);
title('Error de estimación en x');
grid on;
subplot(2,2,3);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),xhat(4,:));
title('Velocidad en y');
grid on;
legend('Real','Estimada');
subplot(2,2,4);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,3)-xhat(4,:)'));
ylim([0 5]);
title('Error de estimación en y');
grid on;

figure(3); %gráfico de la orientación
subplot(1,2,1)
plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),atan2(xhat(7,:),xhat(5,:))/pi*180);
grid on;
title('Orientación');
legend('Real','Estimada');
subplot(1,2,2)
plot(Theta(1:100:end-1,1),mod(360,abs(Theta(1:100:end-1,2)-atan2(xhat(7,:),xhat(5,:))'/pi*180)));
ylim([0 200]);
title('Error de estimación en grados')
grid on;


%% Kalman con sesgo en la velocidad.
% Se tiene el mismo problema que antes, pero ahora la medición tiene un
% sesgo en la velocidad del orden de 1m/s (en cada eje).

clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar-sesgo-vel.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
h = 1; % Período de muestreo.
Abx = Acel(1:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(1:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(1:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz
N = length(Abx); % Cantidad de muestras

% La dimensión del vector $x$ a estimar es n = dim(Pos)+dim(Vel)+dim(Cb), 
% la del ruido $\xi$ del proceso es q = dim(Pos)+dim(Vel)+dim(Cb), y la 
% del ruido $\eta$ de medición es p = dim(Pos)+dim(Vel).
n = size(Preal,2)-1 + size(Vreal,2)-1 + 4;
q = size(Preal,2)-1 + size(Vreal,2)-1 + 4;
p = size(Preal,2)-1 + size(Vreal,2)-1;

% Matrices útiles:
I = eye(size(Preal,2)-1);
O = zeros(size(Preal,2)-1);

% Ahora se tiene un vector z[k], formado por el vector de variables de 
% estado del ejercicio anterior, x[k], y por s[k], el sesgo de la medición.
%
%   z[k+1] = A_moño[k] * z[k] + B_moño[k] * xi[k]
%   y[k]   = C_moño[k] * z[k] + eta[k]
%
%   A_moño[k] = [ A[k] O ]     B_moño = [ B[k] ]     Ck_moño = [ C[k] I ]
%               [ O    I ]              [ O    ]
%        
% Como el sesgo es en velocidad únicamente, la dimensión de s[k] es
% dim(Velocidad) = 2. 

% Dinámica del proceso:
A = zeros(N,n+size(Vradar,2)-1,n+size(Vradar,2)-1);
for i=1:N
    Ak = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
          0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
          0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
          0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
          0 0 0 0 h             w(i)*h         0             0;
          0 0 0 0 -w(i)*h       h              0             0;
          0 0 0 0 0             0              h             w(i)*h;
          0 0 0 0 0             0              -w(i)*h       h];
      
    A(i,:,:) = [Ak                           zeros(length(Ak),length(I));
                zeros(length(I),length(Ak))  I];        
end

% Ruido del proceso:
Bk = [I*h^3/6 O       O       O;
      O       I*h^2/2 O       O;
      O       O       I*h^2/2 O;
      O       O       O       I*h^2/2];
  
B = [Bk;
     O O O O];
        
s_proceso = .01; % desvío del ruido del proceso.
Q = eye(q)*s_proceso;

% Medición. Se mide posición y velocidad con el radar y se guarda cada
% medición en una columna del vector y:
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
C = [I O O O O;  
     O I O O I]; 
  
% Ruido de medición:
sp = 10;
sv = .1;

R = [sp^2*I O;
     O      sv^2*I];
     
% Condiciones iniciales:

% Se elige la estimación inicial de manera que resulte insesgado.
sesgo_vel = 1; 
x0 = [zeros(n,1);
      ones(size(Vradar,2)-1,1)*sesgo_vel];

deltaP = 100;
deltaV = .2;
deltaTheta = 2/9*pi;
deltaC = max(sin(deltaTheta),cos(deltaTheta));
deltaS = sesgo_vel/100;
P = [deltaP*I O        O        O        O;
     O        deltaV*I O        O        O;
     O        O        deltaC*I O        O;
     O        O        O        deltaC*I O;
     O        O        O        O        deltaS*I];
 
% Estimación:
xhat = kalman1(A,B,C,0,0,0,Q,R,x0,P,y); % Estimación por ecuaciones de Riccati
%xhat = kalman2(A,B,C,Q,R,x0,P,y); % Estimación por factorización QR

% Gráficos:
figure(1); %gráfico de la posición
subplot(2,2,1);
plot(Preal(1:100:end-1,1),Preal(1:100:end-1,2),Preal(1:100:end-1,1),xhat(1,:));
title('Posición en x');
grid on;
legend('Real','Estimada');
subplot(2,2,2);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,2)-xhat(1,:)'));
ylim([0 500]);
title('Error de estimación en x');
grid on;
subplot(2,2,3);
plot(Preal(1:100:end-1,1),Preal(1:100:end-1,3),Preal(1:100:end-1,1),xhat(2,:));
title('Posición en y');
grid on;
legend('Real','Estimada');
subplot(2,2,4);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,3)-xhat(2,:)'));
ylim([0 500]);
title('Error de estimación en y');
grid on;


figure(2); %gráfico de la velocidad en x y en y en función del tiempo
subplot(2,2,1);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),xhat(3,:));
title('Velocidad en x');
grid on;
legend('Real','Estimada');
subplot(2,2,2);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,2)-xhat(3,:)'));
ylim([0 5]);
title('Error de estimación en x');
grid on;
subplot(2,2,3);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),xhat(4,:));
title('Velocidad en y');
grid on;
legend('Real','Estimada');
subplot(2,2,4);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,3)-xhat(4,:)'));
ylim([0 5]);
title('Error de estimación en y');
grid on;

figure(3); %gráfico de la orientación
subplot(1,2,1)
plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),atan2(xhat(7,:),xhat(5,:))/pi*180);
grid on;
title('Orientación');
legend('Real','Estimada');
subplot(1,2,2)
plot(Theta(1:100:end-1,1),mod(360,abs(Theta(1:100:end-1,2)-atan2(xhat(7,:),xhat(5,:))'/pi*180)));
ylim([0 200]);
title('Error de estimación en grados')
grid on;


%% Kalman con sesgo en la posición.
% Se tiene el mismo problema que antes, pero ahora la medición tiene un
% sesgo en la velocidad del orden de 1m/s (en cada eje).

clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar-sesgo-pos.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
h = 1; % Período de muestreo.
Abx = Acel(1:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(1:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(1:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz
N = length(Abx); % Cantidad de muestras

% La dimensión del vector $x$ a estimar es n = dim(Pos)+dim(Vel)+dim(Cb), 
% la del ruido $\xi$ del proceso es q = dim(Pos)+dim(Vel)+dim(Cb), y la 
% del ruido $\eta$ de medición es p = dim(Pos)+dim(Vel).
n = size(Preal,2)-1 + size(Vreal,2)-1 + 4;
q = size(Preal,2)-1 + size(Vreal,2)-1 + 4;
p = size(Preal,2)-1 + size(Vreal,2)-1;

% Matrices útiles:
I = eye(size(Preal,2)-1);
O = zeros(size(Preal,2)-1);

% Ahora se tiene un vector z[k], formado por el vector de variables de 
% estado del ejercicio anterior, x[k], y por s[k], el sesgo de la medición.
%
%   z[k+1] = A_moño[k] * z[k] + B_moño[k] * xi[k]
%   y[k]   = C_moño[k] * z[k] + eta[k]
%
%   A_moño[k] = [ A[k] O ]     B_moño = [ B[k] ]     Ck_moño = [ C[k] I ]
%               [ O    I ]              [ O    ]
%        
% Como el sesgo es en velocidad únicamente, la dimensión de s[k] es
% dim(Vposicion) = 2. 

% Dinámica del proceso:
A = zeros(N,n+size(Vradar,2)-1,n+size(Vradar,2)-1);
for i=1:N
    Ak = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
          0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
          0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
          0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
          0 0 0 0 h             w(i)*h         0             0;
          0 0 0 0 -w(i)*h       h              0             0;
          0 0 0 0 0             0              h             w(i)*h;
          0 0 0 0 0             0              -w(i)*h       h];
      
    A(i,:,:) = [Ak                           zeros(length(Ak),length(I));
                zeros(length(I),length(Ak))  I];        
end

% Ruido del proceso:
Bk = [I*h^3/6 O       O       O;
      O       I*h^2/2 O       O;
      O       O       I*h^2/2 O;
      O       O       O       I*h^2/2];
  
B = [Bk;
     O O O O];
        
s_proceso = .01; % desvío del ruido del proceso.
Q = eye(q)*s_proceso;

% Medición. Se mide posición y velocidad con el radar y se guarda cada
% medición en una columna del vector y:
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
C = [I O O O I;  
     O I O O O]; 
  
% Ruido de medición:
sp = 10;
sv = .1;

R = [sp^2*I O;
     O      sv^2*I];
     
% Condiciones iniciales:

% Se elige la estimación inicial de manera que resulte insesgado.
sesgo_pos = 20; 
x0 = [zeros(n,1);
      ones(size(Vradar,2)-1,1)*sesgo_pos];

deltaP = 100;
deltaV = .2;
deltaTheta = 2/9*pi;
deltaC = max(sin(deltaTheta),cos(deltaTheta));
deltaS = sesgo_pos/1000000;
P = [deltaP*I O        O        O        O;
     O        deltaV*I O        O        O;
     O        O        deltaC*I O        O;
     O        O        O        deltaC*I O;
     O        O        O        O        deltaS*I];
 
% Estimación:
xhat = kalman1(A,B,C,0,0,0,Q,R,x0,P,y); % Estimación por ecuaciones de Riccati
%xhat = kalman2(A,B,C,Q,R,x0,P,y); % Estimación por factorización QR

% Gráficos:
figure(1); %gráfico de la posición
subplot(2,2,[1,3]);
plot(Preal(:,2),Preal(:,3),xhat(1,:),xhat(2,:))
title('Posición Estimada');
legend('Real','Estimada');
grid on;
subplot(2,2,2);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,2)-xhat(1,:)'));
ylim([0 100]);
title('Sesgo en x');
grid on;
subplot(2,2,4);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,3)-xhat(2,:)'));
ylim([0 100]);
title('Sesgo en y');
grid on;


figure(2); %gráfico de la velocidad en x y en y en función del tiempo
subplot(2,2,1);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),xhat(3,:));
title('Velocidad en x');
grid on;
legend('Real','Estimada');
subplot(2,2,2);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,2)-xhat(3,:)'));
ylim([0 5]);
title('Error de estimación en x');
grid on;
subplot(2,2,3);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),xhat(4,:));
title('Velocidad en y');
grid on;
legend('Real','Estimada');
subplot(2,2,4);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,3)-xhat(4,:)'));
ylim([0 5]);
title('Error de estimación en y');
grid on;

figure(3); %gráfico de la orientación
subplot(1,2,1)
plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),atan2(xhat(7,:),xhat(5,:))/pi*180);
grid on;
title('Orientación');
legend('Real','Estimada');
subplot(1,2,2)
plot(Theta(1:100:end-1,1),mod(360,abs(Theta(1:100:end-1,2)-atan2(xhat(7,:),xhat(5,:))'/pi*180)));
ylim([0 200]);
title('Error de estimación en grados')
grid on;


%% Kalman extendido

clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
h = 1; % Período de muestreo.
Abx = Acel(100:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(100:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(100:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz
N = length(Abx); % Cantidad de muestras

% La dimensión del vector $x$ a estimar es n = dim(Pos)+dim(Vel)+dim(theta), 
% la del ruido $\xi$ del proceso es p = dim(Pos)+dim(Vel)+dim(theta), y la 
% del ruido $\eta$ de medición es q = dim(Pos)+dim(Vel).
% x = [Px Py Vx Vy theta]'
n = 5;
p = 5;
q = 4;

% Matrices útiles:
I = eye(2); % I[2x2]
O = zeros(2); % O[2x2]

% Dinámica del proceso:

f = @(h,Abx,Aby,w,x) (...
    [x(1)+x(3)*h+(Abx*cos(x(5))-Aby*sin(x(5)))*h^2/2;
     x(2)+x(4)*h+(Abx*cos(x(5))-Aby*sin(x(5)))*h^2/2;
     x(3)+(Abx*cos(x(5))-Aby*sin(x(5)))*h;
     x(4)+(Abx*sin(x(5))-Aby*cos(x(5)))*h;
     x(5)+w*h]);

A = @(h,Abx,Aby,w,xhat) (...
    [1 0 h 0 -(Abx*sin(xhat(5))+Aby*cos(xhat(5)))*h^2/2;
     0 1 0 h (Abx*cos(xhat(5))-Aby*sin(xhat(5)))*h^2/2; 
     0 0 1 0 -(Abx*sin(xhat(5))+Aby*cos(xhat(5)))*h;
     0 0 0 1 (Abx*cos(xhat(5))-Aby*sin(xhat(5)))*h;
     0 0 0 0 1]);

% Ruido del proceso:
B = [h^3/6   0     0     0     0;
     0       h^3/6 0     0     0;
     0       0     h^2/2 0     0;
     0       0     0     h^2/2 0;
     0       0     0     0     h^2/2];
 
s_proceso = .00001; % desvío del ruido del proceso.
Q = eye(p)*s_proceso;

% Medición. Se mide posición y velocidad con el radar y se guarda cada
% medición en una columna del vector y:
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
C = [I O zeros(2,1);
     O I zeros(2,1)]; 

% Ruido de medición:
sp = 10;
sv = .1;

R = [sp^2*I O;   
     O      sv^2*I];

% Condiciones iniciales:
deltaP = 100;
deltaV = .2;
deltaTheta = 2/9*pi;
%deltaC = max(sin(deltaTheta),cos(deltaTheta));
P = [deltaP*I   O          [0 0]';
     O          deltaV*I   [0 0]';
     [0 0]      [0 0]      deltaTheta];
x0 = zeros(n,1);

xhat = zeros(n,size(y,2));
xhat(:,1) = x0;

for k=1:N-1
    % Predicción
    Ak = A(h,Abx(k),Aby(k),w(k),xhat(:,k));
    xhat(:,k+1) = f(h,Abx(k),Aby(k),w(k),xhat(:,k));
    P = Ak*P*Ak' + B*Q*B';
    % Actualización:
    K = P*C'*(inv(R+C*P*C'));
    xhat(:,k+1) = xhat(:,k+1) + K*(y(:,k+1)-C*xhat(:,k+1));
    P = (eye(n) - K*C)*P;    
end

% Gráficos:
figure(1); %gráfico de la posición
subplot(2,2,1);
plot(Preal(1:100:end-1,1),Preal(1:100:end-1,2),Preal(1:100:end-1,1),xhat(1,:));
title('Posición en x');
grid on;
legend('Real','Estimada');
subplot(2,2,2);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,2)-xhat(1,:)'));
ylim([0 500]);
title('Error de estimación en x');
grid on;
subplot(2,2,3);
plot(Preal(1:100:end-1,1),Preal(1:100:end-1,3),Preal(1:100:end-1,1),xhat(2,:));
title('Posición en y');
grid on;
legend('Real','Estimada');
subplot(2,2,4);
plot(Preal(1:100:end-1,1),abs(Preal(1:100:end-1,3)-xhat(2,:)'));
ylim([0 500]);
title('Error de estimación en y');
grid on;


figure(2); %gráfico de la velocidad en x y en y en función del tiempo
subplot(2,2,1);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),xhat(3,:));
title('Velocidad en x');
grid on;
legend('Real','Estimada');
subplot(2,2,2);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,2)-xhat(3,:)'));
ylim([0 5]);
title('Error de estimación en x');
grid on;
subplot(2,2,3);
plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),xhat(4,:));
title('Velocidad en y');
grid on;
legend('Real','Estimada');
subplot(2,2,4);
plot(Vreal(1:100:end-1,1),abs(Vreal(1:100:end-1,3)-xhat(4,:)'));
ylim([0 5]);
title('Error de estimación en y');
grid on;

figure(3); %gráfico de la orientación
subplot(1,2,1)
plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),xhat(5,:)/pi*180+360);
grid on;
title('Orientación');
legend('Real','Estimada');
subplot(1,2,2)
plot(Theta(1:100:end-1,1),mod(360,abs(Theta(1:100:end-1,2)-xhat(5,:)'/pi*180)));
ylim([0 200]);
title('Error de estimación en grados')
grid on;



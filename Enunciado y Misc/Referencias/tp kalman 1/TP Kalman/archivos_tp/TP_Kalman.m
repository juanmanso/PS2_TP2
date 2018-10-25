%% FILTRO DE KALMAN
%
% Con el objetivo de usar el filtro de Kalman que relaciona las variables
% de estado en forma lineal, se consideró un primer modelo del sistema como
% el siguiente: las variables de estado son la posición $P^e$ del vehículo con respecto
% a las coordenadas terrestres $x$, $y$, la velocidad $V^e$ del vehículo
% con respecto a las mismas coordenadas y el conjunto de variables:
%
% $$
% \begin{align}
% c_{11} &= \cos(\theta)
% c_{21} &= -\sin(\theta)
% c_{12} &= \sin(\theta)
% c_{22} &= \cos(\theta)
% \end{align}
% $$
%
% donde $\theta$ es la orientación del vehículo. De esta forma, la dinámica
% del sistema puede describirse mediante el siguiente conjunto de
% ecuaciones:
% 
% (poner ecuaciones)

%%
tic
clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
Ts = 1;
Abx = Acel(1:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(1:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(1:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz


% A = @(Abx,Aby,w,h) ([1 0 h 0 Abx/2*h^2  Aby/2*h^2   0          0;
%                      0 1 0 h 0          0           Abx/2*h^2  Aby/2*h^2;
%                      0 0 1 0 Abx*h      Aby*h       0          0;
%                      0 0 0 1 0          0           Abx*h      Aby*h;
%                      0 0 0 0 h          w*h         0          0;
%                      0 0 0 0 -w*h       h           0          0;
%                      0 0 0 0 0          0           h          w*h;
%                      0 0 0 0 0          0           -w*h       h]);

A = zeros(length(Abx),8,8);
h = Ts;
for i=1:length(Abx)
    A(i,:,:) = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
                0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
                0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
                0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
                0 0 0 0 h             w(i)*h         0             0;
                0 0 0 0 -w(i)*h       h              0             0;
                0 0 0 0 0             0              h             w(i)*h;
                0 0 0 0 0             0              -w(i)*h       h];
end
    

% Matrices para ecuaciones del sistema:
I = eye(2);
O = zeros(2);
% A = f(tk) = depende del tiempo (se define en las ecuaciones de Riccati)

B = [I O O O;
     O I O O;
     O O I O;
     O O O I]; 



C = [I O O O;
     O I O O]; % Se mide posición y velocidad con el radar.
 
sp = 10;
sv = .1;
R = [sp*I O;   
     O    sv*I];
 
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
% Se grafica la posición, velocidad (en x y) y la orientación del vehículo.
n = length(Abx)-1; % cantidad de muestras del proceso
for j=0:3
    Q = [I*h/6 O     O   O;
         O     I*h/2 O   O;
         O     O     I*h O;
         O     O     O   I*h]*10^(-j); % Ruido del proceso
    % Condiciones iniciales:
    %xhat = zeros(length(B),n);
    deltaP = 100;
    deltaV = .2;
    deltaTheta = 2/9*pi;
    deltaC = max(sin(deltaTheta),cos(deltaTheta));
    P = [deltaP*I O        O        O;
         O        deltaV*I O        O;
         O        O        deltaC*I O;
         O        O        O        deltaC*I];
    h = Ts;

    xhat = kalman1(A,B,C,Q,R,zeros(length(P),1),P,y);
    figure(1); %gráfico de la posición
    subplot(2,2,j+1); hold on;
    plot(Preal(:,2),Preal(:,3),xhat(1,:),xhat(2,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(2); %gráfico de la velocidad en x en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),xhat(3,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(3); %gráfico de la velocidad en y en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),xhat(4,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(4); %gráfico de la orientación
    subplot(2,2,j+1); hold on;
    plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),atan2(xhat(7,:),xhat(5,:))/pi*180);
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
end
toc

%% Ejercicio 2:
% Se mide la posición y la velocidad con el radar, y se quiere estimar la 
% posición, velocidad y orientación del vehículo. Para verificar la
% estimación se tienen Preal, Vreal y Theta.
tic
clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
Ts = 1;
Abx = Acel(1:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(1:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(1:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz


% Matrices para ecuaciones del sistema:
I = eye(2);
O = zeros(2);
% A = f(tk) = depende del tiempo (se define en las ecuaciones de Riccati)

B = [I O O O;
     O I O O;
     O O I O;
     O O O I]; 

h = Ts/1000;

C = [I O O O;
     O I O O]; % Se mide posición y velocidad con el radar.
 
sp = 10;
sv = .1;
R = [sp*I O;   
     O    sv*I];
 
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
% Se grafica la posición, velocidad (en x y) y la orientación del vehículo.
n = length(Abx)-1; % cantidad de muestras del proceso
for j=0:3
    Q = [I*h/6 O     O   O;
         O     I*h/2 O   O;
         O     O     I*h O;
         O     O     O   I*h]*10^(-j); % Ruido del proceso
    % Condiciones iniciales:
    xhat = zeros(length(B),n);
    deltaP = 100;
    deltaV = .2;
    deltaTheta = 2/9*pi;
    deltaC = max(sin(deltaTheta),cos(deltaTheta));
    P = [deltaP*I O        O        O;
         O        deltaV*I O        O;
         O        O        deltaC*I O;
         O        O        O        deltaC*I];
    h = Ts;
    
    for i=1:n-1            
        A = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
             0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
             0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
             0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
             0 0 0 0 h             w(i)*h         0             0;
             0 0 0 0 -w(i)*h       h              0             0;
             0 0 0 0 0             0              h             w(i)*h;
             0 0 0 0 0             0              -w(i)*h       h];
        xhat(:,i+1) = A*xhat(:,i); % Predicción
        P = A*P*A' + B*Q*B';       
        K = P*C'*(inv(R+C*P*C'));
        xhat(:,i+1) = xhat(:,i+1) + K*(y(:,i+1) -C*xhat(:,i+1)); % Actualización
        P = (eye(length(A)) - K*C)*P;
    end
    figure(1); %gráfico de la posición
    subplot(2,2,j+1); hold on;
    plot(Preal(:,2),Preal(:,3),xhat(1,:),xhat(2,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(2); %gráfico de la velocidad en x en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),xhat(3,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(3); %gráfico de la velocidad en y en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),xhat(4,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(4); %gráfico de la orientación
    subplot(2,2,j+1); hold on;
    plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),atan2(xhat(7,:),xhat(5,:))/pi*180);
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
end
toc

%% Ejercicio 3:
% El problema se mantiene igual, pero ahora el radar tiene un sesgo en la
% velocidad de aproximadamente 1m/s en cada eje.

clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar-sesgo-vel.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
Ts = 1;
Abx = Acel(1:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(1:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(1:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz


% Matrices para ecuaciones del sistema:
% Ahora se tiene un vector z = [x; s] con s el sesgo en velocidad.
I = eye(2);
O = zeros(2);
% A = f(tk) = depende del tiempo (se define en las ecuaciones de Riccati)

B = [I O O O;
     O I O O;
     O O I O;
     O O O I;
     O O O O;
     O O O I]; 

h = Ts/1000;

C = [I O O O I O;
     O I O O O I]; % Se mide posición y velocidad con el radar con un sesgo en velocidad:
                   % C_moño = [C I]
 
sp = 10;
sv = .1;
R = [sp*I O;   
     O    sv*I];
 
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
% Se grafica la posición, velocidad (en x y) y la orientación del vehículo.
n = length(Abx)-1; % cantidad de muestras del proceso

for j=0:3
    Q = [I*h/6 O     O   O;
         O     I*h/2 O   O;
         O     O     I*h O;
         O     O     O   I*h]*10^(-j); % Ruido del proceso
    % Condiciones iniciales:
    zhat = zeros(length(B),n);
    zhat(9:12,1) = [0 0 1 1]';
    deltaP = 100;
    deltaV = .2;
    deltaTheta = 2/9*pi;
    deltaC = max(sin(deltaTheta),cos(deltaTheta));
    P = [deltaP*I O        O        O        O O;
         O        deltaV*I O        O        O O;
         O        O        deltaC*I O        O O;
         O        O        O        deltaC*I O O;
         O        O        O        O        O O;
         O        O        O        O        O O];
    h = Ts;
    
    for i=1:n-1            
        A1 = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
              0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
              0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
              0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
              0 0 0 0 h             w(i)*h         0             0;
              0 0 0 0 -w(i)*h       h              0             0;
              0 0 0 0 0             0              h             w(i)*h;
              0 0 0 0 0             0              -w(i)*h       h];
        A = [A1                  zeros(length(A1),4);
             zeros(4,length(A1)) eye(4)               ];
        zhat(:,i+1) = A*zhat(:,i); % Predicción
        P = A*P*A' + B*Q*B';       
        K = P*C'*(inv(R+C*P*C'));
        zhat(:,i+1) = zhat(:,i+1) + K*(y(:,i+1) -C*zhat(:,i+1)); % Actualización
        P = (eye(length(A)) - K*C)*P;
    end
    figure(1); %gráfico de la posición
    subplot(2,2,j+1); hold on;
    plot(Preal(:,2),Preal(:,3),zhat(1,:),zhat(2,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(2); %gráfico de la velocidad en x en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),zhat(3,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(3); %gráfico de la velocidad en y en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),zhat(4,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(4); %gráfico de la orientación
    subplot(2,2,j+1); hold on;
    plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),atan2(zhat(7,:),zhat(5,:))/pi*180);
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
end

%% Ejercicio 5: Kalman por factorización de Cholesky y QR

clear;
close all;

load('Acel.mat'); % Me devuelve Acel, muestreada a 100Hz
load('Gyro.mat'); % Me devuelve Gyro, muestreada a 100Hz
load('Radar.mat'); % Me devuelve Pradar y Vradar, muestreadas a 1Hz
load('trayectoria.mat'); % Me devuelve Preal, Vreal y Theta, muestreadas a 100Hz

% Se submuestrea la aceleración y la velocidad angular a 1Hz.
Ts = 1;
Abx = Acel(1:100:end,2); % Aceleración en x, en coordenadas no inerciales, submuestreada a 1Hz
Aby = Acel(1:100:end,3); % Aceleración en y, en coordenadas no inerciales, submuestreada a 1Hz
w = Gyro(1:100:end,2); % Velocidad angular en z, en coordenadas no inerciales, submuestreada a 1Hz


% Matrices para ecuaciones del sistema:
I = eye(2);
O = zeros(2);
h = Ts/1000;
% A = f(tk) = depende del tiempo (se define en las ecuaciones de Riccati)

B = [I O O O;
     O I O O;
     O O I O;
     O O O I]; 

C = [I O O O;
     O I O O]; % Se mide posición y velocidad con el radar.
 
sp = 10;
sv = .1;
R = [sp*I O;   
     O    sv*I];
 
y = [Pradar(:,2)'; 
     Pradar(:,3)';
     Vradar(:,2)';
     Vradar(:,3)'];
 
% Se grafica la posición, velocidad (en x y) y la orientación del vehículo.
N = length(Abx)-1; % cantidad de muestras del proceso
Q = [I*h/6 O     O   O;
     O     I*h/2 O   O;
     O     O     I*h O;
     O     O     O   I*h]*10^(-1); % Ruido del proceso

% Condiciones iniciales:
xhat = zeros(length(B),N);
deltaP = 100;
deltaV = .2;
deltaTheta = 2/9*pi;
deltaC = max(sin(deltaTheta),cos(deltaTheta));
P = [deltaP*I O        O        O;
     O        deltaV*I O        O;
     O        O        deltaC*I O;
     O        O        O        deltaC*I];
h = Ts;
    
% Predicción inicial:
A = [1 0 h 0 Abx(1)/2*h^2  Aby(1)/2*h^2   0             0;
     0 1 0 h 0             0              Abx(1)/2*h^2  Aby(1)/2*h^2;
     0 0 1 0 Abx(1)*h      Aby(1)*h       0             0;
     0 0 0 1 0             0              Abx(1)*h      Aby(1)*h;
     0 0 0 0 h             w(1)*h         0             0;
     0 0 0 0 -w(1)*h       h              0             0;
     0 0 0 0 0             0              h             w(1)*h;
     0 0 0 0 0             0              -w(1)*h       h];
xhat(:,1) = A*xhat(:,1);
P = A*P*A' + B*Q*B'; 
P_sqrt = chol(P);
R_sqrt = chol(R);
Q_sqrt = chol(Q);
M = [R_sqrt                     C*P_sqrt zeros(length(R),length(Q));
     zeros(length(Q),length(R)) A*P_sqrt B*Q_sqrt                 ];
[q, r] = qr(M');
r
Z = r(length(R)+1:end,length(R)+1:end)';

fQ = [-y(:,2)'*inv(R_sqrt') xhat(:,1)'*inv(P_sqrt') zeros(1,length(Q))]*q;
W = fQ(length(R)+1:length(R)+length(A)+1);
size(Z)
size(W)
xhat(:,2) = Z*W';
%%
for i=2:N-1
    M = [R_sqrt                     C*Z zeros(length(R),length(Q));
         zeros(length(Q),length(R)) A*Z B*Q_sqrt                 ];
    [q, r] = qr(M');
    Z = r(lenght(R)+1:end,length(R)+1:end)';
    fQ = [-y(:,i)'*inv(R_sqrt') W zeros(1,length(Q))]*q;
    W = fQ(length(R)+1:length(R)+length(A)+1);
    xhat(i+1) = Z*W';
end

figure(1);
plot(Preal(:,2),Preal(:,3),xhat(1,:),xhat(2,:))
grid on;
legend('Real','Estimada');

%%     




    for i=2:N-1
        A = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
             0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
             0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
             0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
             0 0 0 0 h             w(i)*h         0             0;
             0 0 0 0 -w(i)*h       h              0             0;
             0 0 0 0 0             0              h             w(i)*h;
             0 0 0 0 0             0              -w(i)*h       h];
         R_sqrt = chol(R);
         Q_sqrt = chol(Q);
         M = [R_sqrt                     C*P_sqrt zeros(length(R),length(Q));
              zeros(length(Q),length(R)) A*P_sqrt B*Q_sqrt                 ];
         [q, r] = qr(M');
         Z = r(length(A)+1:length(A)+1+length(P),length(A)+1:end)';
         %f = [-y(:,i)'*inv(R_sqrt') xhat(:,i)*inv(P_sqrt') zeros(length(R),length(Q))]*q;
         %size(-y(:,i)'*inv(R_sqrt'))
         %size(W2)
         %size(zeros(1,length(Q)))
         f = [-y(:,i)'*inv(R_sqrt') W2 zeros(1,length(Q))]*q;
         %f(:,length(A)+1:length(A)+1+length(Z))
         W2 = f(:,length(A)+1:length(A)+1+length(Z));
         size(Z)
         size(W2)
         xhat(:,i+1) = Z*W2';
    end
         
             
    
%     for i=1:n-1            
%         A = [1 0 h 0 Abx(i)/2*h^2  Aby(i)/2*h^2   0             0;
%              0 1 0 h 0             0              Abx(i)/2*h^2  Aby(i)/2*h^2;
%              0 0 1 0 Abx(i)*h      Aby(i)*h       0             0;
%              0 0 0 1 0             0              Abx(i)*h      Aby(i)*h;
%              0 0 0 0 h             w(i)*h         0             0;
%              0 0 0 0 -w(i)*h       h              0             0;
%              0 0 0 0 0             0              h             w(i)*h;
%              0 0 0 0 0             0              -w(i)*h       h];
%         xhat(:,i+1) = A*xhat(:,i); % Predicción
%         P = A*P*A' + B*Q*B';       
%         K = P*C'*(inv(R+C*P*C'));
%         xhat(:,i+1) = xhat(:,i+1) + K*(y(:,i+1) -C*xhat(:,i+1)); % Actualización
%         P = (eye(length(A)) - K*C)*P;
%     end
    figure(1); %gráfico de la posición
    subplot(2,2,j+1); hold on;
    plot(Preal(:,2),Preal(:,3),xhat(1,:),xhat(2,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(2); %gráfico de la velocidad en x en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,2),Vreal(1:100:end-1,1),xhat(3,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(3); %gráfico de la velocidad en y en función del tiempo
    subplot(2,2,j+1); hold on;
    plot(Vreal(1:100:end-1,1),Vreal(1:100:end-1,3),Vreal(1:100:end-1,1),xhat(4,:))
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
    figure(4); %gráfico de la orientación
    subplot(2,2,j+1); hold on;
    plot(Theta(:,1),Theta(:,2),Theta(1:100:end-1,1),atan2(xhat(7,:),xhat(5,:))/pi*180);
    title(sprintf('Q = I/%d',10^j),'Fontsize',9);
    grid on;
    legend('Real','Estimada');
    hold off;
end








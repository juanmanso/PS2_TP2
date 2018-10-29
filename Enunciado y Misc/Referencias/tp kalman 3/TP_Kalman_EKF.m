clear all
close all
clc

load './archivos_tp/Acel.mat'
load './archivos_tp/Gyro.mat'
load './archivos_tp/Radar.mat'
load './archivos_tp/trayectoria.mat'

% Filtro de Kalman extendido

N = length(Acel);
Ax_b = Acel(:,2);
Ay_b = Acel(:,2);
w = Gyro(:,2);

T_k = 0.01; % Paso 10ms

I = eye(2);
O = 0*I;

A = zeros(5,5,N);
 
B = eye(5);

C = [[I O [0;0]];
	 [O I [0;0]]];

sigma_zeta = 0.1; % Hago la prueba con un valor chico
Q = (sigma_zeta)^2*eye(5); % Q arbitraria

sigma_pos = 100; % +-100m
sigma_vel = 0.1; % +-0.1m/s
R = [[(sigma_pos)^2*I  O];
     [O               (sigma_vel)^2*I]];
 
% Inicializacion
X = zeros(5,N);
x = [0 0 0 0 0]'; % x0 arbitrario

% 40° = 0.6981 rad
sigma_tita = 0.6981;
     
P = [[(sigma_pos)^2*I O                     [0;0]];
     [O               (sigma_vel)^2*I       [0;0]];
     [0  0             0   0              (sigma_tita)^2]];

 
 
for i = 1:1:N
    A(:,:,i) = [[1 0 T_k 0 0];
                [0 1 0 T_k 0];
                [0 0 1 0   -Ax_b(i)*sin(x(5))-Ay_b(i)*cos(x(5))];
                [0 0 0 1    Ax_b(i)*cos(x(5))-Ay_b(i)*sin(x(5))];
                [0 0 0 0        1]];

    % Prediccion
    x = [[x(1)+T_k*x(3)];
         [x(2)+T_k*x(4)];
         [x(3)+T_k*(Ax_b(i)*cos(x(5))-Ay_b(i)*sin(x(5)))];
         [x(4)+T_k*(Ax_b(i)*sin(x(5))+Ay_b(i)*cos(x(5)))];
         [x(5)+T_k*w(i)]];
    P = A(:,:,i)*P*A(:,:,i)' + B*Q*B';
    
    if(mod(i,100) == 0)
        % Actualizacion
        K = P*C'*((R + C*P*C')^(-1));
        x = x + K*...
  ([Pradar(i/100,2);Pradar(i/100,3);Vradar(i/100,2);Vradar(i/100,3)]-C*x);
        P = (eye(5) - K*C)*P;
    end
    
    X(:,i) = x;
end

% Posicion
figure(1)
plot(X(1,:),X(2,:))
hold on
plot(Preal(:,2),Preal(:,3),'r')

% Velocidad
figure(2)
plot(X(3,:),X(4,:))
hold on
plot(Vreal(:,2),Vreal(:,3),'r')

% Orientacion
figure(3)
plot(rad2deg(X(5,:)))
hold on
plot(Theta(:,2),'r') % Esta en grados de -180 a 180

% Error en posicion
figure(4)
plot((X(1,:)'-Preal(:,2)).^2+(X(2,:)'-Preal(:,3)).^2)

% Error en velocidad
figure(5)
plot((X(3,:)'-Vreal(:,2)).^2+(X(4,:)'-Vreal(:,3)).^2)

% Error en orientacion
figure(6)
plot(rad2deg(X(5,:))'-Theta(:,2))

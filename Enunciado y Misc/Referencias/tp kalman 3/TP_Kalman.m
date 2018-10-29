clear all
close all
clc

load './archivos_tp/Acel.mat'
load './archivos_tp/Gyro.mat'
load './archivos_tp/Radar.mat'
load './archivos_tp/trayectoria.mat'

% Filtro de Kalman, extendiendo el vector de estados

N = length(Acel);
Ax_b = Acel(:,2);
Ay_b = Acel(:,2);
w = Gyro(:,2);

T_k = 0.01; % Paso 10ms

I = eye(2);
O = 0*I;

A = zeros(6,6,N);

for i = 1:1:N
A(:,:,i) = [[1 0 T_k 0 (T_k^2)/2*Ax_b(i) (T_k^2)/2*Ay_b(i)];
            [0 1 0 T_k (T_k^2)/2*Ay_b(i) -(T_k^2)/2*Ax_b(i)];
            [0 0 1 0    T_k*Ax_b(i)        T_k*Ax_b(i)];
            [0 0 0 1    T_k*Ay_b(i)      -T_k*Ax_b(i)];
            [0 0 0 0        1               T_k*w(i)];
            [0 0 0 0        -T_k*w(i)           1]];
end
 
B = eye(6);

C = [[I O O];
	 [O I O]];

sigma_zeta = 0.1; % Hago la prueba con un valor chico
Q = (sigma_zeta)^2*eye(6); % Q arbitraria

sigma_pos = 100; % +-100m
sigma_vel = 0.1; % +-0.1m/s
R = [[(sigma_pos)^2*I  O];
     [O               (sigma_vel)^2*I]];
 
% Inicializacion
X = zeros(6,N);
x = [0 0 0 0 1 0]'; % x0 arbitrario

% 40° = 0.6981 rad
sigma_tita = 0.6981;
P_aux = [[(cos(sigma_tita))^2 0];
         [0               (sin(sigma_tita))^2]];
     
P = [[(sigma_pos)^2*I O                     O];
     [O               (sigma_vel)^2*I       O];
     [O               O                 P_aux]];

 
for i = 1:1:N
    % Prediccion
    x = A(:,:,i)*x;
    P = A(:,:,i)*P*A(:,:,i)' + B*Q*B';
    
    if(mod(i,100) == 0)
        % Actualizacion
        K = P*C'*((R + C*P*C')^(-1));
        x = x + K*...
  ([Pradar(i/100,2);Pradar(i/100,3);Vradar(i/100,2);Vradar(i/100,3)]-C*x);
        P = (eye(6) - K*C)*P;
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
plot(atan2d(-X(6,:),X(5,:)))
hold on
plot(Theta(:,2),'r') % Esta en grados de -180 a 180

% Posicion y orientacion
figure(4)
x = 0:.05:2*pi;
y = sin(x);
z = zeros(size(x));
col = atan2d(-X(6,:),X(5,:));% This is the color, vary with x in this case.
surface([X(1,:);X(1,:)],[X(2,:);X(2,:)],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);

% Posicion y velocidad    
figure(5)
%quiver(X(1,1:100),X(2,1:100),X(3,1:100),X(4,1:100))
quiver(Preal(:,2),Preal(:,3),Vreal(:,2),Vreal(:,3),0.01)


% Error en posicion
figure(6)
plot((X(1,:)'-Preal(:,2)).^2+(X(2,:)'-Preal(:,3)).^2)

% Error en velocidad
figure(7)
plot((X(3,:)'-Vreal(:,2)).^2+(X(4,:)'-Vreal(:,3)).^2)

% Error en orientacion
figure(8)
plot(atan2d(-X(6,:),X(5,:))'-Theta(:,2))

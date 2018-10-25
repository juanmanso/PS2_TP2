clear all;
close all;

%parametros

dP_i = 100;
dV_i = 0.2;
dTheta_g_i = 40;
dTheta_i = dTheta_g_i*pi/180;
dP_rad = 10;
dV_rad = 0.1;

load('Acel.mat');
load('Gyro.mat');
load('Radar.mat');
load('trayectoria.mat');

T_s_D = 0.01;
T_s_S = 1;
T_s_rat = round(T_s_S/T_s_D);

%modelo tomando 0 como variable de estado (sistema no lineal)

% estado : [Px Py Vx Vy O Asx Asy]
syms A_x A_y w_t O A_x_s A_y_s

aux_Vx_O = -sin(O) * (A_x) - cos(O) * (A_y);
aux_Vy_O =  cos(O) * (A_x) - sin(O) * (A_y);

aux_mat = [ aux_Vx_O;
            aux_Vy_O;   ];

A_t = [     eye(2), T_s_D * eye(2), zeros(2,1);
            zeros(2), eye(2), T_s_D * aux_mat;
            zeros(1,2), zeros(1,2), eye(1);];

B_t = 0;

C_t = [ eye(2), zeros(2), zeros(2,1);
        zeros(2), eye(2), zeros(2,1);  ];
    
%kalman extendido

theta_i = Theta(1,2)*pi/180 + normrnd(0,dTheta_i);

X_0 = [ Preal(1,2)+ normrnd(0,dP_i);
        Preal(1,3)+ normrnd(0,dP_i);
        Vreal(1,2)+ normrnd(0,dV_i);
        Vreal(1,3)+ normrnd(0,dV_i);
        theta_i;];

    
P_0 = diag( [ dP_i^2, dP_i^2, dV_i^2, dV_i^2, (dTheta_i)^2] );

Q = 0;
R = diag( [dP_rad^2, dP_rad^2, dV_rad^2, dV_rad^2] );

aux = size(Acel);
samples_size = aux(1);

aux = size(X_0);
X_k = zeros(aux(1), aux(2), samples_size);
X_kp1_k = zeros(aux(1), aux(2), samples_size);
aux = size(P_0);
P_k = zeros(aux(1), aux(2), samples_size);
P_kp1_k = zeros(aux(1), aux(2), samples_size);
aux = size(A_t);
A_k = zeros(aux(1), aux(2), samples_size);
aux = size(B_t);
B_k = zeros(aux(1), aux(2), samples_size);
aux = size(C_t);
C_k = zeros(aux(1), aux(2), samples_size);
aux = size(C_t * X_0);
G_k = zeros(aux(1), aux(2), samples_size);

for i = 1:samples_size
    
    if i == 1
        X_est = X_0;
        P_est = P_0;
    else
        X_est = X_k(:,:,i-1);
        P_est = P_k(:,:,i-1);
    end

    A_x = Acel(i,2);
    A_y = Acel(i,3);
    w_t = Gyro(i,2);
    P_x = X_est(1);
    P_y = X_est(2);
    V_x = X_est(3);
    V_y = X_est(4);
    O = X_est(5);

    aux_Vx_O = -sin(O) * (A_x) - cos(O) * (A_y);
    aux_Vy_O =  cos(O) * (A_x) - sin(O) * (A_y);

    aux_mat = [ aux_Vx_O;
                aux_Vy_O;   ];

    A = [     eye(2), T_s_D * eye(2), zeros(2,1);
              zeros(2), eye(2), T_s_D * aux_mat;
              zeros(1,2), zeros(1,2), eye(1);];

    B = 0;

    C = [ eye(2), zeros(2), zeros(2,1);
          zeros(2), eye(2), zeros(2,1);  ];
    
    X_pre = [   P_x + T_s_D * V_x;
                P_y + T_s_D * V_y;
                V_x + T_s_D * ((A_x) * cos(O) - (A_y) * sin(O));
                V_y + T_s_D * ((A_x) * sin(O) + (A_y) * cos(O));
                O + T_s_D * w_t; ];
    
    P_pre = A * P_est * A' + B * Q * B';
    
    X_k(:,:,i) = X_pre;
    P_k(:,:,i) = P_pre;
    
    if rem(i,T_s_rat) == 0
        
        j = i/T_s_rat;
        Y = [   Pradar(j,2);
                Pradar(j,3);
                Vradar(j,2);
                Vradar(j,3);   ];
        K = P_pre * C' * inv(R + C * P_pre * C');
        G =  Y - C * X_pre;
        X_est = X_pre + K * G;
        P_est = (eye(size(P_pre)) - K * C) * P_pre;
        
        X_k(:,:,i) = X_est;
        P_k(:,:,i) = P_est;
        G_k(:,:,j) = G;        
    end
        
    A_k(:,:,i)= A;
    B_k(:,:,i)= B;
    C_k(:,:,i)= C;
    
end

figure
hold on;
plot(Preal(:,2),Preal(:,3))
plot(Pradar(:,2),Pradar(:,3))
plot(reshape(X_k(1,:,:),1,[]),reshape(X_k(2,:,:),1,[]))
legend('real','radar','estimado')

 waitforbuttonpress
 
figure
hold on;
plot3(Vreal(:,2),Vreal(:,3),1:samples_size)
plot3(Vradar(:,2),Vradar(:,3),100:100:samples_size)
plot3(reshape(X_k(3,:,:),1,[]),reshape(X_k(4,:,:),1,[]),1:samples_size)
legend('real','radar','estimado')

 waitforbuttonpress

%comparo estados estimados con reales

 figure
hold on;
 plot(1:samples_size,Preal(:,2))
 plot(1:samples_size,reshape(X_k(1,:,:),[],1))
 
 waitforbuttonpress

 figure
hold on;
 plot(1:samples_size,Preal(:,3))
 plot(1:samples_size,reshape(X_k(2,:,:),[],1))
 
 waitforbuttonpress

 figure
hold on;
 plot(1:samples_size,Vreal(:,2))
 plot(1:samples_size,reshape(X_k(3,:,:),[],1))
 
 waitforbuttonpress

 figure
hold on;
 plot(1:samples_size,Vreal(:,3))
 plot(1:samples_size,reshape(X_k(4,:,:),[],1))

 waitforbuttonpress
 
figure
hold on;
 plot(1:samples_size,Theta(:,2)*pi/180)
 plot(1:samples_size,reshape(X_k(5,:,:),[],1))

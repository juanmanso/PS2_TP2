%parametros

dP_i = 100;
dV_i = 0.2;
dTheta_g_i = 40;
dTheta_i = dTheta_g_i*pi/180;
dVs_i = 1;
dP_rad = 10;
dV_rad = 0.1;


load('Radar.mat');
PradarReal = Pradar;
VradarReal = Vradar;
load('Acel.mat');
load('gyro.mat');
load('Radar-sesgo-vel.mat');
load('trayectoria.mat');
Vsesgo = Vradar - VradarReal;

T_s_D = 0.01;
T_s_S = 1;
T_s_rat = round(T_s_S/T_s_D);

%modelo tomando cos(0) y sin(0) como variable de estado (sistema lineal)

% estado : [Px Py Vx Vy cos(0) sin(0) Vsx Vsy]
syms A_x_t A_y_t w_t

A_b_t =  [  A_x_t;
            A_y_t;  ];

aux_A = [   A_x_t, -A_y_t;
            A_y_t, A_x_t;   ];
  
aux_w = [   0, -w_t;
            w_t, 0;   ]; 

A_con_t = [ zeros(2), eye(2), zeros(2), zeros(2);
            zeros(2), zeros(2), aux_A, zeros(2);
            zeros(2), zeros(2), aux_w, zeros(2);
            zeros(2), zeros(2), zeros(2), zeros(2);];
        
A_disc_t_aprox = eye(size(A_con_t)) + T_s_D * A_con_t;
A_disc_t_exacta = expm(A_con_t * T_s_D);

A_t = A_disc_t_aprox;

B_t = 0;

C_t = [ eye(2), zeros(2), zeros(2), zeros(2);
        zeros(2), eye(2), zeros(2), eye(2);   ];
    
% kalman simple

theta_i = Theta(1,2)*pi/180 + normrnd(0,dTheta_i);

X_0 = [ Preal(1,2)+ normrnd(0,dP_i);
        Preal(1,3)+ normrnd(0,dP_i);
        Vreal(1,2)+ normrnd(0,dV_i);
        Vreal(1,3)+ normrnd(0,dV_i);
        cos(theta_i);
        sin(theta_i);
        dVs_i;
        dVs_i;];
    
P_0 = diag( [ dP_i^2, dP_i^2, dV_i^2, dV_i^2, (dTheta_i)^4, (dTheta_i)^2, dVs_i^2, dVs_i^2] );

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
    
%     esto es lo correcto pero matlab es muy lento asi que harcodeo matrices
%     A = subs(A_t, [ A_x_t, A_y_t, w_t ], [Acel(i,2), Acel(i,3) , Gyro(i,2)]);
%     B = subs(B_t, [ A_x_t, A_y_t, w_t ], [Acel(i,2), Acel(i,3) , Gyro(i,2)]); 
%     C = subs(C_t, [ A_x_t, A_y_t, w_t ], [Acel(i,2), Acel(i,3) , Gyro(i,2)]); 

    A_x_t = Acel(i,2);
    A_y_t = Acel(i,3);
    w_t = Gyro(i,2);
 
% discretizacion aproximada (primer orden)
%
    A = [   [ 1, 0, 1/100,     0,         0,          0, 0, 0]
            [ 0, 1,     0, 1/100,         0,          0, 0, 0]
            [ 0, 0,     1,     0, A_x_t/100, -A_y_t/100, 0, 0]
            [ 0, 0,     0,     1, A_y_t/100,  A_x_t/100, 0, 0]
            [ 0, 0,     0,     0,         1,   -w_t/100, 0, 0]
            [ 0, 0,     0,     0,   w_t/100,          1, 0, 0]
            [ 0, 0,     0,     0,         0,          0, 1, 0]
            [ 0, 0,     0,     0,         0,          0, 0, 1]  ];
        
    B = 0; 
    
    C = [   [ 1, 0, 0, 0, 0, 0, 0, 0]
            [ 0, 1, 0, 0, 0, 0, 0, 0]
            [ 0, 0, 1, 0, 0, 0, 1, 0]
            [ 0, 0, 0, 1, 0, 0, 0, 1] ];
    
    if i == 1
        X_est = X_0;
        P_est = P_0;
    else
        X_est = X_k(:,:,i-1);
        P_est = P_k(:,:,i-1);
    end
    
    X_pre = A * X_est;
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

figure
hold on;
plot3(Vreal(:,2),Vreal(:,3),1:samples_size)
plot3(Vradar(:,2),Vradar(:,3),100:100:samples_size)
plot3(reshape(X_k(3,:,:),1,[]),reshape(X_k(4,:,:),1,[]),1:samples_size)
legend('real','radar','estimado')

%comparo estados estimados con reales
% 
%  figure
% hold on;
%  plot(1:samples_size,Preal(:,2))
%  plot(1:samples_size,reshape(X_k(1,:,:),[],1))
%  
%  waitforbuttonpress
% 
%  figure
% hold on;
%  plot(1:samples_size,Preal(:,3))
%  plot(1:samples_size,reshape(X_k(2,:,:),[],1))
%  
%  waitforbuttonpress
% 
%  figure
% hold on;
%  plot(1:samples_size,Vreal(:,2))
%  plot(1:samples_size,reshape(X_k(3,:,:),[],1))
%  
%  waitforbuttonpress
% 
%  figure
% hold on;
%  plot(1:samples_size,Vreal(:,3))
%  plot(1:samples_size,reshape(X_k(4,:,:),[],1))
% 
%  waitforbuttonpress
%  
% figure
% hold on;
%  plot(1:samples_size,cos(Theta(:,2)*pi/180))
%  plot(1:samples_size,reshape(X_k(5,:,:),[],1))
%  
%  waitforbuttonpress
% 
%  figure
% hold on;
%  plot(1:samples_size,sin(Theta(:,2)*pi/180))
%  plot(1:samples_size,reshape(X_k(6,:,:),[],1))
% 
%  
%  waitforbuttonpress
% 
%  figure
% hold on;
%  plot(100:100:samples_size,Vsesgo(:,2))
%  plot(1:samples_size,reshape(X_k(7,:,:),[],1))
% 
%  waitforbuttonpress
% 
%  figure
% hold on;
%  plot(100:100:samples_size,Vsesgo(:,3))
%  plot(1:samples_size,reshape(X_k(8,:,:),[],1))



config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ5 - Kalman Con QR / Cholesky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acel_str = load('Acel.mat');
gyro_str = load('Gyro.mat');
radar_str = load('Radar.mat');
radar_b_v_str = load('Radar-sesgo-vel.mat');
radar_b_p_str = load('Radar-sesgo-pos.mat');
tray_str = load('trayectoria.mat');

% Definición de los datos.
Acel = acel_str.Acel;
Gyro = gyro_str.Gyro;
Pradar = radar_str.Pradar;
Vradar = radar_str.Vradar;
Pradar_pb = radar_b_p_str.Pradar;
Vradar_pb = radar_b_p_str.Vradar;
Pradar_vb = radar_b_v_str.Pradar;
Vradar_vb = radar_b_v_str.Vradar;
Preal = tray_str.Preal;
Vreal = tray_str.Vreal;
Theta = tray_str.Theta;

% Definición de frecuencias y períodos de muestreo.
Fs_gyro_acel=100;
Ts_gyro_acel=1/Fs_gyro_acel;
Ts=Ts_gyro_acel;

% Cantidad de mediciones y estados.
dim=2;
cant_mediciones=dim*2;
cant_estados=dim*4;
cant_muestras=length(Acel);
cant_muestras_rad=length(Pradar);

% Varianzas.
var_xip=2*(0.01^3)/6;  % CHEQUEAR QUÉ VALORES PONER ACÁ
var_xiv=0.01^2;
var_xic=0.01^2;
sigma_etap=10;
sigma_etav=0.1;

% Definiciones de matrices y vectores de Kalman.

% Matriz de rotación.
Cbe = [cos(Theta(2:2:end,2)) -sin(Theta(2:2:end,2)); sin(Theta(2:2:end,2)) cos(Theta(2:2:end,2))];

I = eye(dim); 
O = zeros(dim);

Bk1 = eye(cant_estados);

C = [I O O O;
     O I O O];

Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,2*dim)*var_xic]);

R = diag([ones(1,dim)*sigma_etap^2 ones(1,dim)*sigma_etav^2]);

yk = [Pradar(:,2:3) Vradar(:,2:3)];
yk = kron(yk, [1;zeros(100 - 1, 1)]);

% Inicio del Kalman.
cov_p = [1 1]*100;
cov_v = [1 1]*0.2;
cov_c = [1 1]*1;

% Condiciones iniciales.
x0 = [100 100 0.2 0.2 cos(40*pi/180) -sin(40*pi/180) sin(40*pi/180) cos(40*pi/180)]';
P0_0 = diag([cov_p, cov_v, cov_c, cov_c]);

x = x0;
P = P0_0;

xk1_k1 = x;
Pk1_k1 = P;
g = yk(1,:)';
    
for i = 1 : cant_muestras - 1
    
    %A Discreta
    Ad= [1 0 Ts 0  (Acel(i,1)*Ts^2)/2 (Acel(i,2)*Ts^2)/2 0                   0;
        0 1 0  Ts 0                  0                  (Acel(i,1)*Ts^2)/2 (Acel(i,2)*Ts^2)/2;
        0 0 1  0  Acel(i,1)*Ts       Acel(i,2)*Ts       0                   0;
        0 0 0  1  0                  0                  Acel(i,1)*Ts        Acel(i,2)*Ts;
        0 0 0  0  Ts                 Gyro(i,2)*Ts       0                   0;
        0 0 0  0  -Gyro(i,2)*Ts      Ts                 0                   0;
        0 0 0  0  0                  0                  Ts                  Gyro(i,2)*Ts;
        0 0 0  0  0                  0                  -Gyro(i,2)*Ts       Ts];

    %Algoritmo QR
	[xk1_k1, Pk1_k1, gk] = kalman_qr(xk1_k1, Pk1_k1, Ad, Bk1, C, Qd, R, yk(i,:));
    
	g = [g gk];
    x = [x xk1_k1];

end

%Grafico de medida, estimada, ruidosa
h1=figure;
hold on;
grid;
plot(Preal(:,2),Preal(:,3),'r','LineWidth',2);
plot(x(1,:),x(2,:),'--g','LineWidth',2);
title('Estimación de la trayectoria');
if(EsMatlab == 1)
    legend('Real','Estimada','location','SouthEast');
    xlabel('Posición x');
    ylabel('Posición y');
else
    legend(['Real';'Estimada'],'location','SouthEast');
    xlabel('Posicion $x$ [\si{\m}]');
    ylabel('Posicion $y$ [\si{\m}]');
end

	% Grafico del estado posición en función del tiempo
	h2=figure;
    subplot(3,1,[1 2]);
    hold on;
	grid;
    plot(Preal(:,2),'LineWidth',2);
	plot(x(1,:),'--','LineWidth',2);
    plot(Preal(:,3),'LineWidth',2);
	plot(x(2,:),'--','LineWidth',2);
    title('Estados de posición');
    xlim([0 length(x(1,:))]);
if(EsMatlab == 1)
    legend('X Real','X Estimada','Y Real','Y Estimada','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posición');
else
    legend(['X Real';'X Estimada';'Y Real';'Y Estimada'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posición');
    end
    subplot(3,1,3);
    hold on;
    grid;
    plot(Preal(:,2)-x(1,:)','LineWidth',1);
    plot(Preal(:,3)-x(2,:)','LineWidth',1);
    xlim([0 length(x(1,:))]);
    if(EsMatlab == 1)
    legend('Error X','Error Y','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posición');
else
    legend(['Error X';'Error Y'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posición');
    end
    h2.Position=[0 0 700 1200];
    h2.PaperUnits='points';
    h2.PaperSize=[700 1200];
	
% Grafico del estado velocidad en función del tiempo
	h2=figure;
	hold on;
	grid;
    plot(Vreal(:,2),'LineWidth',2);
	plot(x(3,:),'--','LineWidth',2);
    plot(Vreal(:,3),'LineWidth',2);
	plot(x(4,:),'--','LineWidth',2);
	title('Estados de velocidad');
    if(EsMatlab == 1)
    legend('X Real','X Estimada','Y Real','Y Estimada','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
else
    legend(['X Real';'X Estimada';'Y Real';'Y Estimada'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
    end

    
% Grafico de tita función del tiempo
	h3=figure;
	hold on;
	grid;
    plot(Theta(:,2),'LineWidth',2);
	plot((180/pi)*acos(x(5,:)).*(-sign(x(6,:))),'--','LineWidth',2);
	title('Estados de theta');
    if(EsMatlab == 1)
    legend('\theta Real','\theta Estimada','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
else
    legend(['\theta Real';'\theta Estimada'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
end

	% GrÃ¡fico de correlaciÃ³n de innovaciones (debe ser ruido blanco)
	covx_g = xcorr(g(1,:)');
	covy_g = xcorr(g(2,:)');
	
	figure
	plot(covx_g)
	grid
	title('Covarianza innovaciones x')
	
	figure
	plot(covy_g)
	grid
	title('Covarianza innovaciones y')
	
	% Observabilidad
	Obs = obsv(Ad,C);
	rango_obs = rank(Obs);
	estados_no_observables = cant_estados - rango_obs
	
	

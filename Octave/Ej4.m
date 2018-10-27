config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TP 2 KALMAN EJ 4 - Estimación con sesgo en posición a partir de mediciones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acel_str = load('archivos_tp\Acel.mat');
gyro_str = load('archivos_tp\Gyro.mat');
radar_str = load('archivos_tp\Radar.mat');
radar_b_v_str = load('archivos_tp\Radar-sesgo-vel.mat');
radar_b_p_str = load('archivos_tp\Radar-sesgo-pos.mat');
tray_str = load('archivos_tp\trayectoria.mat');

% Definición de los datos
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

% Definición de frecuencias y períodos de muestreo
Fs_gyro_acel=100;
Ts_gyro_acel=1/Fs_gyro_acel;
Ts=Ts_gyro_acel;

% Cantidad de mediciones y estados
dim=2;
cant_mediciones=dim*2;
cant_estados=dim*4;
cant_muestras=length(Acel);
cant_muestras_rad=length(Pradar);

% Varianzas
var_xip=100;  % CHEQUEAR QUÉ VALORES PONER ACÁ
var_xiv=10;
var_xic=1;
sigma_etap=10;
sigma_etav=0.1;

% Definiciones de matrices y vectores de Kalman
    % Matriz de rotaciÃ³n
    Cbe = [cos(Theta(2:2:end,2)) -sin(Theta(2:2:end,2)); sin(Theta(2:2:end,2)) cos(Theta(2:2:end,2))];

    I=eye(dim); 
    O=zeros(dim);

    % Variable de estado
    % X=[P;V;C11;C12;C21;C22], C11=Cbe(1,1), C12=Cbe(1,2), C21=Cbe(2,1), C22=Cbe(2,2)

    %A= [0 0 1 0 0 0 0 0;
    %    0 0 0 1 0 0 0 0;
    %    0 0 0 0 ax ay 0 0;
    %    0 0 0 0 0 0 ax ay;
    %    0 0 0 0 0 wb 0 0;
    %    0 0 0 0 -wb 0 0 0;
    %    0 0 0 0 0 0 0 wb;
    %    0 0 0 0 0 0 -wb 0];
   
    Bk1=eye(cant_estados);
    C=[I O O O;
       O I O O];
    M_eta = [randn(dim,cant_muestras)*10;
            randn(dim,cant_muestras)*0.1;
            zeros(dim,cant_muestras);
            zeros(dim,cant_muestras)];
    Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,2*dim)*var_xic]);
    R = diag([ones(1,dim)*sigma_etap^2 ones(1,dim)*sigma_etav^2]);
    yk=[Pradar_pb(:,2:3)+randn(cant_muestras_rad,dim)*sigma_etap Vradar_pb(:,2:3)+randn(cant_muestras_rad,dim)*sigma_etav];
    
    cov_p = [1 1]*100;
	cov_v = [1 1]*0.2;
	cov_c = [1 1]*1;
	
	x0 = [100 100 0.2 0.2 cos(40*pi/180) -sin(40*pi/180) sin(40*pi/180) cos(40*pi/180)]';
	P0_0 = diag([cov_p, cov_v, cov_c, cov_c]);

	x = x0;
	P = P0_0;
	xk1_k1 = x;
	Pk1_k1 = P;
	g = yk(1,:)';
for i=1:cant_muestras_rad
%A discreta
Ad= [1 0 Ts 0  (Acel(i,1)*Ts^2)/2 (Acel(i,2)*Ts^2)/2 0                   0;
    0 1 0  Ts 0                  0                  (Acel(i,1)*Ts^2)/2 (Acel(i,2)*Ts^2)/2;
    0 0 1  0  Acel(i,1)*Ts       Acel(i,2)*Ts       0                   0;
    0 0 0  1  0                  0                  Acel(i,1)*Ts        Acel(i,2)*Ts;
    0 0 0  0  Ts                 Gyro(i,2)*Ts       0                   0;
    0 0 0  0  -Gyro(i,2)*Ts      Ts                 0                   0;
    0 0 0  0  0                  0                  Ts                  Gyro(i,2)*Ts;
    0 0 0  0  0                  0                  -Gyro(i,2)*Ts       Ts];

% PredicciÃ³n
		xk_k1 = Ad * xk1_k1;
		Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
		gk = [innovaciones(yk(i,:),C,xk_k1)];
%	
%		% CorrecciÃ³n
		Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
		xk_k = xk_k1 + Kk*(gk);
		Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
%		
%		% ActualizaciÃ³n
		xk1_k1 = xk_k;
		Pk1_k1 = Pk_k;
%	
%	
%		% Guardo
		g = [g gk];
		x = [x xk_k];
		P = [P; Pk_k];

end

%Grafico de medida, estimada, ruidosa
figure;
hold on;
grid;
plot(Preal(:,2),Preal(:,3),'LineWidth',3);
plot(x(1,:),x(2,:),'r','LineWidth',2);
title('Estimación de la trayectoria');
legend('Real','Estimada','location','SouthEast');
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
	figure;
	hold on;
	grid;
    plot(Preal(:,2),'LineWidth',2);
	plot(1:100:cant_muestras,x(1,:),'--','LineWidth',2);
    plot(Preal(:,3),'LineWidth',2);
	plot(1:100:cant_muestras,x(2,:),'--','LineWidth',2);
	title('Estados de posición');
    if(EsMatlab == 1)
    legend('X Real','X Estimada','Y Real','Y Estimada','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posición');
else
    legend(['X Real';'X Estimada';'Y Real';'Y Estimada'],'location','SouthEast');
    xlabel('Posicion $x$ [\si{\m}]');
    ylabel('Posicion $y$ [\si{\m}]');
    end
	
% Grafico del estado velocidad en función del tiempo
	figure;
	hold on;
	grid;
    plot(Vreal(:,2),'LineWidth',2);
	plot(1:100:cant_muestras,x(3,:),'--','LineWidth',2);
    plot(Vreal(:,3),'LineWidth',2);
	plot(1:100:cant_muestras,x(4,:),'--','LineWidth',2);
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
	figure;
	hold on;
	grid;
    plot(Theta(:,2),'LineWidth',2);
	plot(1:100:cant_muestras,(180/pi)*acos(x(5,:)).*(-sign(x(6,:))),'--','LineWidth',2);
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
    
    

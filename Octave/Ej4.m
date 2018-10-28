config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TP 2 KALMAN EJ 4 - Estimaci�n con sesgo en posici�n a partir de mediciones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acel_str = load('archivos_tp\Acel.mat');
gyro_str = load('archivos_tp\Gyro.mat');
radar_str = load('archivos_tp\Radar.mat');
radar_b_v_str = load('archivos_tp\Radar-sesgo-vel.mat');
radar_b_p_str = load('archivos_tp\Radar-sesgo-pos.mat');
tray_str = load('archivos_tp\trayectoria.mat');

% Definici�n de los datos
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

% Definici�n de frecuencias y per�odos de muestreo
Fs_gyro_acel=100;
Ts_gyro_acel=1/Fs_gyro_acel;
Ts=Ts_gyro_acel;

% Cantidad de mediciones y estados
dim=2;
cant_mediciones=dim*2;
cant_estados=dim*4+2;  % Se a�aden dos por el c�lculo del sesgo!
cant_muestras=length(Acel);
cant_muestras_rad=length(Pradar);

% Varianzas
var_xip=0.001;  % CHEQUEAR QU� VALORES PONER AC�
var_xiv=0.1;
var_xic=0.01;
var_xib=0.01^10;
sigma_etap=10;
sigma_etav=0.1;

% Definiciones de matrices y vectores de Kalman
    % Matriz de rotación
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
   
    Bk1=diag([ones(1,dim), ones(1,dim), ones(1,dim),ones(1,dim), zeros(1,dim)]);
    C=[I O O O I;
       O I O O O];
    Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,2*dim)*var_xic,ones(1,dim)*var_xib]);
    R = diag([ones(1,dim)*sigma_etap^2 ones(1,dim)*sigma_etav^2]);
    yk=[Pradar_pb(:,2:3) Vradar_pb(:,2:3)];
    yk=kron(yk,[1;zeros(100-1,1)]);
    
    cov_p = [1 1]*100;
	cov_v = [1 1]*0.2;
	cov_c = [1 1]*1;
	cov_b = [1 1]*50;
    
	x0 = [100 100 0.2 0.2 cos(40*pi/180) -sin(40*pi/180) sin(40*pi/180) cos(40*pi/180) 20 20]';
	P0_0 = diag([cov_p, cov_v, cov_c, cov_c, cov_b]);

	x = x0;
	P = P0_0;
	xk1_k1 = x;
	Pk1_k1 = P;
	g = yk(1,:)';
for i=1:cant_muestras-1
%A discreta

    ad35=Acel(i,2)*Ts-(Acel(i,3)*Gyro(i,2)*Ts^2)/2;
    ad36=Acel(i,3)*Ts+(Acel(i,2)*Gyro(i,2)*Ts^2)/2;
    ad47=Acel(i,2)*Ts-(Acel(i,3)*Gyro(i,2)*Ts^2)/2;
    ad48=Acel(i,3)*Ts+(Acel(i,2)*Gyro(i,2)*Ts^2)/2;

    Ad= [1 0 Ts 0 (Acel(i,2)*Ts^2)/2     (Acel(i,3)*Ts^2)/2     0                      0;
        0 1 0  Ts 0                      0                      (Acel(i,2)*Ts^2)/2    (Acel(i,3)*Ts^2)/2;
        0 0 1  0  ad35                   ad36                   0                      0;
        0 0 0  1  0                      0                      ad47                   ad48;
        0 0 0  0  1-((Gyro(i,2)*Ts)^2)/2 Gyro(i,2)*Ts           0                      0;
        0 0 0  0  -Gyro(i,2)*Ts          1-((Gyro(i,2)*Ts)^2)/2 0                      0;
        0 0 0  0  0                      0                      1-((Gyro(i,2)*Ts)^2)/2 Gyro(i,2)*Ts;
        0 0 0  0  0                      0                      -Gyro(i,2)*Ts          1-((Gyro(i,2)*Ts)^2)/2];

Ad=[Ad, zeros(size(Ad,1),dim);
    zeros(dim, size(Ad,2)),  I];

% Predicción
		xk_k1 = Ad * xk1_k1;
		Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
% 		gk = [innovaciones(yk(i,:),C,xk_k1)];
%	
%		% Correcci�n
if(yk(i,1)~=0)
        gk = [innovaciones(yk(i,:),C,xk_k1)];
		Kk = Pk_k1 * C'*(R+ C*Pk_k1*C')^-1;
		xk_k = xk_k1 + Kk*(gk);
		Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
%		
%		% Actualizaci�n
		xk1_k1 = xk_k;
		Pk1_k1 = Pk_k;
        
else
        % Actualizaci�n
		xk1_k1 = xk_k1;
		Pk1_k1 = Pk_k1;
end
%	
%	
%		% Guardo
		g = [g gk];
        x = [x xk1_k1]; %Guardo el que acabo de actualizar para la siguiente iteraci�n
% 		x = [x xk_k];
% 		P = [P; Pk_k]; % Desactivado para no llenar memoria, no se utiliza

end

%Grafico de medida, estimada, ruidosa
h1=figure;
hold on;
grid;
plot(Preal(:,2),Preal(:,3),'r','LineWidth',2);
plot(x(1,:),x(2,:),'--g','LineWidth',2);
title('Estimaci�n de la trayectoria');
if(EsMatlab == 1)
    legend('Real','Estimada','location','SouthEast');
    xlabel('Posici�n x');
    ylabel('Posici�n y');
else
    legend(['Real';'Estimada'],'location','SouthEast');
    xlabel('Posicion $x$ [\si{\m}]');
    ylabel('Posicion $y$ [\si{\m}]');
end

	% Grafico del estado posici�n en funci�n del tiempo
	h2=figure;
    subplot(3,1,[1 2]);
    hold on;
	grid;
    plot(Preal(:,2),'LineWidth',2);
	plot(x(1,:),'--','LineWidth',2);
    plot(Preal(:,3),'LineWidth',2);
	plot(x(2,:),'--','LineWidth',2);
    title('Estados de posici�n');
    xlim([0 length(x(1,:))]);
if(EsMatlab == 1)
    legend('X Real','X Estimada','Y Real','Y Estimada','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posici�n');
else
    legend(['X Real';'X Estimada';'Y Real';'Y Estimada'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posici�n');
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
    ylabel('Posici�n');
else
    legend(['Error X';'Error Y'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posici�n');
    end
    h2.Position=[0 0 700 1200];
    h2.PaperUnits='points';
    h2.PaperSize=[700 1200];
	
% Grafico del estado velocidad en funci�n del tiempo
	h3=figure;
    subplot(3,1,[1 2]);
	hold on;
	grid;
    plot(Vreal(:,2),'LineWidth',2);
	plot(x(3,:),'--','LineWidth',2);
    plot(Vreal(:,3),'LineWidth',2);
	plot(x(4,:),'--','LineWidth',2);
	title('Estados de velocidad');
    xlim([0 length(x(3,:))]);
    if(EsMatlab == 1)
    legend('X Real','X Estimada','Y Real','Y Estimada','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
else
    legend(['X Real';'X Estimada';'Y Real';'Y Estimada'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
    end
    subplot(3,1,3);
    hold on;
    grid;
    plot(Vreal(:,2)-x(3,:)','LineWidth',1);
    plot(Vreal(:,3)-x(4,:)','LineWidth',1);
    xlim([0 length(x(3,:))]);
    if(EsMatlab == 1)
    legend('Error X','Error Y','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posici�n');
else
    legend(['Error X';'Error Y'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Posici�n');
    end
    h3.Position=[0 0 700 1200];
    h3.PaperUnits='points';
    h3.PaperSize=[700 1200];

    
% Grafico de tita funci�n del tiempo
	h3=figure;
	hold on;
	grid;
    plot(cos(Theta(:,2)*pi/180),'LineWidth',2);
	plot(x(5,:),'--','LineWidth',2);
	title('Estados cos(\theta(t))');
    if(EsMatlab == 1)
    legend('\theta Real','\theta Estimada','location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
else
    legend(['\theta Real';'\theta Estimada'],'location','SouthEast');
    xlabel('Tiempo');
    ylabel('Velocidad');
end

	% Gráfico de correlación de innovaciones (debe ser ruido blanco)
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
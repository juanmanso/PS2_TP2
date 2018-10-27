config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ KALMAN - Estimación a partir de mediciones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acel_str = load('archivos_tp\Acel.mat');
gyro_str = load('archivos_tp\Gyro.mat');
radar_str = load('archivos_tp\Radar.mat');
radar_b_v_str = load('archivos_tp\Radar-sesgo-vel.mat');
radar_b_p_str = load('archivos_tp\Radar-sesgo-pos.mat');
tray_str = load('archivos_tp\trayectoria.mat');

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

Fs_gyro_acel=100;
Ts_gyro_acel=1/Fs_gyro_acel;


% Matriz de rotación
Cbe = [cos(Theta(2:2:end,2)) -sin(Theta(2:2:end,2)); sin(Theta(2:2:end,2)) cos(Theta(2:2:end,2))];

dim=2;
cant_mediciones=dim*4;
I=eye(dim);
O=zeros(dim);

% Cantidad de muestras
cant_muestras=length(Acel);
cant_muestras_rad=length(Pradar);

% Variable de estado
% X=[P;V;C11;C12;C21;C22], C11=Cbe(1,1), C12=Cbe(1,2), C21=Cbe(2,1),
% C22=Cbe(2,2)

%A= [0 0 1 0 0 0 0 0;
%    0 0 0 1 0 0 0 0;
%    0 0 0 0 ax ay 0 0;
%    0 0 0 0 0 0 ax ay;
%    0 0 0 0 0 wb 0 0;
%    0 0 0 0 -wb 0 0 0;
%    0 0 0 0 0 0 0 wb;
%    0 0 0 0 0 0 -wb 0];

Ts=Ts_gyro_acel;
var_xip=100;  % CHEQUEAR QU� VALORES PONER AC�
var_xiv=10;
var_xic=1;
sigma_etap=10;
sigma_etav=0.1;

M_eta = [randn(dim,cant_muestras)*10;
        randn(dim,cant_muestras)*0.1;
        zeros(dim,cant_muestras);
       	zeros(dim,cant_muestras)];
Bk1=eye(8);
C=eye(8);
Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,2*dim)*var_xic]);

yk=[Pradar(:,2:3)+randn(cant_muestras_rad,dim)*sigma_etap Vradar(:,2:3)+randn(cant_muestras_rad,dim)*sigma_etav Gyro(1:100:end-1,:)];
% 	% Datos
% 	var_xip = 3e-4;
% 	var_xiv = 2e-3;
% 	var_xia = 1e-2;

    cov_p = [1 1]*100;
	cov_v = [1 1]*0.2;
	cov_c = [1 1]*1;
	
	x0 = [0 0 0 0 0 0 0 0]';
	P0_0 = diag([cov_p, cov_v, cov_c, cov_c]);

	x = x0;
	P = P0_0;
	xk1_k1 = x;
	Pk1_k1 = P;
	g = yk(1,:)';
for i=1:cant_muestras
%A discreta
Ad= [1 0 Ts 0  (Acel(i,1)*Ts^2)/2 (Acel(i,2)*Ts^2)/2 0                   0;
    0 1 0  Ts 0                  0                  (Acel(i,1)*Ts^2)/2 (Acel(i,2)*Ts^2)/2;
    0 0 1  0  Acel(i,1)*Ts       Acel(i,2)*Ts       0                   0;
    0 0 0  1  0                  0                  Acel(i,1)*Ts        Acel(i,2)*Ts;
    0 0 0  0  Ts                 Gyro(i,2)*Ts       0                   0;
    0 0 0  0  -Gyro(i,2)*Ts      Ts                 0                   0;
    0 0 0  0  0                  0                  Ts                  Gyro(i,2)*Ts;
    0 0 0  0  0                  0                  -Gyro(i,2)*Ts       Ts];

% Predicción
		xk_k1 = Ad * xk1_k1;
		Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
		gk = [innovaciones(yk(i,:),C,xk_k1)];
%	
%		% Corrección
		Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
		xk_k = xk_k1 + Kk*(gk);
		Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
%		
%		% Actualización
		xk1_k1 = xk_k;
		Pk1_k1 = Pk_k;
%	
%	
%		% Guardo
		g = [g gk];
		x = [x xk_k];
		P = [P; Pk_k];

end
return;
%	
%	dim = 2;			% Se considera sólo x e y
%	tipos_variables = 3;		% Posición, Velocidad, Aceleración
%	cant_mediciones = length(Pos);
%	cant_estados = tipos_variables * dim;
%	
%	bool_p = 1;
%	bool_v = 0;
%	bool_a = 0;
%	
%	%%%%%%%%%%%%%%
%	%%% 1a Defina las variables de estado
%	%%%%%%%%%%%%%%%
%	
%	%%%%%%%%%%%%%%
%	%%% 1b Ad y Q_d
%	%%%%%%%%%%%%%%%
%	
%	% Datos
%	var_xip = 3e-4;
%	var_xiv = 2e-3;
%	var_xia = 1e-2;
%	
%	%%%
%	T = Tiempo(2:end)-Tiempo(1:end-1);
%	T = 1;
%	
%	% Variable de estado X = [P;V;A]
%	I = eye(dim);
%	Ad =	[I	I.*T	(T.^2)/2.*I;
%		 I*0	I	T.*I;
%		 I*0	I*0	I;];
%	
%	Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,dim)*var_xia]); %Sólo para x e y
%	
%	
%	
%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	% EJ 2
%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	
%	cov_p = [1 1]*100^2;
%	cov_v = [1 1]*1;
%	cov_a = [1 1]*0.1;
%	
%	x0 = [40 -200 0 0 0 0]';
%	P0_0 = diag([cov_p, cov_v, cov_a]);
%	
%	%% a)
%	%%%%% y_k = [I 0 0] [pk vk ak]' + ruido \eta
%	sigma_etap = 60;
%	sigma_etav = 2;
%	sigma_etaa = 0.1;
%	
%	%%% Para hacer AWGN, randn(fila,col)*sigma_etap
%	
%	Bk1 = eye(cant_estados);
%	
%	% C = [eye(dim)*bool_p eye(dim)*bool_v eye(dim)*bool_a];	% Obsoleto
%	C =	[eye(dim*bool_p) zeros(dim*bool_p) zeros(dim*bool_p);
%		 zeros(dim*bool_v) eye(dim*bool_v) zeros(dim*bool_v);
%		 zeros(dim*bool_a) zeros(dim*bool_a) eye(dim*bool_a)];
%	
%	M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p; 
%		randn(dim,cant_mediciones)*sigma_etav*bool_v;
%	       	randn(dim,cant_mediciones)*sigma_etaa*bool_a];
%	
%	yk = C * [Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)]' + (C*M_eta);
%	yk = yk'; % Así tiene la forma de Pos
%	
%	R = diag([ones(1,dim*bool_p)*sigma_etap^2 ones(1,dim*bool_v)*sigma_etav^2 ones(1,dim*bool_a)*sigma_etaa^2]);
%	
%	
%	%%% ALGORITMO %%%%
%	x = x0;
%	P = P0_0;
%	xk1_k1 = x;
%	Pk1_k1 = P;
%	g = yk(1,:)';
%	
%	for i=1:cant_mediciones-1
%		% Predicción
%		xk_k1 = Ad * xk1_k1;
%		Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
%		gk = [innovaciones(yk(i,:),C,xk_k1)];
%	
%		% Corrección
%		Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
%		xk_k = xk_k1 + Kk*(gk);
%		Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
%		
%		% Actualización
%		xk1_k1 = xk_k;
%		Pk1_k1 = Pk_k;
%	
%	
%		% Guardo
%		g = [g gk];
%		x = [x xk_k];
%		P = [P; Pk_k];
%	end
%	
%	
%	% Grafico de medida, estimada, ruidosa
%	figure
%	hold on
%	grid
%	plot(x(1,:),x(2,:),'LineWidth',3)
%	plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
%	plot(yk(:,1),yk(:,2),'color',myGreen)
%	title('Estimación de la trayectoria - Medición de $p$');
%	legend(['Estimada';'Medida';'Ruidosa'],'location','SouthEast');
%	xlabel('Posición $x$ [\si{\m}]');
%	ylabel('Posición $y$ [\si{\m}]');
%	return;
%	
%	% Grafico del estado posición en función del tiempo
%	figure
%	hold on
%	grid
%	plot(x(1,:),'LineWidth',2)
%	plot(x(2,:),'color',myGreen,'LineWidth',2)
%	title('Estados de posición');
%	
%	
%	% Grafico del estado velocidad en función del tiempo
%	figure
%	hold on
%	grid
%	plot(x(3,:),'LineWidth',2)
%	plot(x(4,:),'color',myGreen,'LineWidth',2)
%	title('Estados de velocidad');
%	
%	
%	% Grafico del estado aceleración en función del tiempo
%	figure
%	hold on
%	grid
%	plot(x(5,:),'LineWidth',2)
%	plot(x(6,:),'color',myGreen,'LineWidth',2)
%	title('Estados de aceleración');
%	
%	
%	% Gráfico de correlación de innovaciones (debe ser ruido blanco)
%	covx_g = xcorr(g(1,:)');
%	covy_g = xcorr(g(2,:)');
%	
%	figure
%	plot(covx_g)
%	grid
%	title('Covarianza innovaciones x')
%	
%	figure
%	plot(covy_g)
%	grid
%	title('Covarianza innovaciones y')
%	
%	% Observabilidad
%	Obs = obsv(Ad,C);
%	rango_obs = rank(Obs);
%	estados_no_observables = cant_estados - rango_obs
%	
%	

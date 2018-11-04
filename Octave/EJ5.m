config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ5 - Kalman Con QR / Cholesky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imprimir imágenes?
bool_print=0;

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
% var_xip=2*(0.01^3)/6;  % CHEQUEAR QUÉ VALORES PONER ACÁ
% var_xiv=0.01^2;
% var_xic=0.01^4;

var_xip=2*(0.01^3)/6;  % CHEQUEAR QUÉ VALORES PONER ACÁ
% var_xip=0.01^5;
var_xiv=0.01^2;
var_xic=0.01^4;

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
yk = kron(yk, [zeros(100 - 1, 1);1]);

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

    %Algoritmo QR
	[xk1_k1, Pk1_k1, gk] = kalman_qr(xk1_k1, Pk1_k1, Ad, Bk1, C, Qd, R, yk(i,:));
    
    if (yk(i,:) ~= 0)
        g = [g gk];
    end
    x = [x xk1_k1];

end

%Grafico de medida, estimada, ruidosa
h1=figure;
hold on;
grid;
plot(Preal(:,2),Preal(:,3),'LineWidth',2,'Color',colors.greenery);
plot(x(1,:),x(2,:),'--r','LineWidth',2);
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

     setGraphSize(h1,wide_1);

if bool_print
    print('../Informe/Figuras/graf_ej5','-dpdf','-bestfit');
end

	% Grafico del estado posición en función del tiempo
	h2=figure;
    subplot(3,1,[1 2]);
    hold on;
	grid;
    plot(0:Ts:300,Preal(:,2),'LineWidth',2,'Color',colors.carrot);
	plot(0:Ts:300,x(1,:),'--','LineWidth',2,'Color',colors.blithe);
    plot(0:Ts:300,Preal(:,3),'LineWidth',2,'Color',colors.raspberry);
	plot(0:Ts:300,x(2,:),'--','LineWidth',2,'Color',colors.castaway);
    title('Estados de posición');
    xlim([0 300]);
if(EsMatlab == 1)
    legend('X Real','X Estimada','Y Real','Y Estimada','location','SouthEast');
    xlabel('Tiempo (s)');
    ylabel('Posición (m)');
else
    legend(['X Real';'X Estimada';'Y Real';'Y Estimada'],'location','SouthEast');
    xlabel('Tiempo (s)');
    ylabel('Posición (m)');
end
    subplot(3,1,3);
    hold on;
    grid;
    plot(0:Ts:300,Preal(:,2)-x(1,:)','LineWidth',1);
    plot(0:Ts:300,Preal(:,3)-x(2,:)','LineWidth',1);
    xlim([0 300]);
    if(EsMatlab == 1)
    legend('Error X','Error Y','location','SouthEast');
    xlabel('Tiempo (s)');
    ylabel('Error (m)');
else
    legend(['Error X';'Error Y'],'location','SouthEast');
    xlabel('Tiempo (s)');
    ylabel('Error (m)');
    end
    setGraphSize(h2,vertical_1);
    
if bool_print
    print('../Informe/Figuras/graf_ej5_pos','-dpdf','-bestfit');
end

	
% Grafico del estado velocidad en función del tiempo
	h3=figure;
	hold on;
	grid;
    plot(0:Ts:300,Vreal(:,2),'LineWidth',2,'Color',colors.carrot);
	plot(0:Ts:300,x(3,:),'--','LineWidth',2,'Color',colors.blithe);
    plot(0:Ts:300,Vreal(:,3),'LineWidth',2,'Color',colors.raspberry);
	plot(0:Ts:300,x(4,:),'--','LineWidth',2,'Color',colors.castaway);
	title('Estados de velocidad');
    xlim([0 300]);
    if(EsMatlab == 1)
    legend('X Real','X Estimada','Y Real','Y Estimada','location','SouthEast');
    xlabel('Tiempo (s)');
    ylabel('Velocidad (m/s)');
else
    legend(['X Real';'X Estimada';'Y Real';'Y Estimada'],'location','SouthEast');
    xlabel('Tiempo (s)');
    ylabel('Velocidad (m/s)');
    end
    setGraphSize(h3,wide_1);

if bool_print
    print('../Informe/Figuras/graf_ej5_vel','-dpdf','-bestfit');
end

% Grafico de Cnn en función del tiempo
	h4=figure;
    subplot 311;
	hold on;
	grid;
    plot(0:Ts:300,cos(Theta(:,2)*pi/180),'LineWidth',3,'Color',colors.mint);
	plot(0:Ts:300,x(5,:),'g','LineWidth',2,'Color',colors.gecko);
    plot(0:Ts:300,x(8,:),'--','LineWidth',2,'Color',colors.raspberry);
    xlim([0 300]);
    if(EsMatlab == 1)
    legend('cos(\theta) Real','C_{11} Estimada','C_{22} Estimada','location','NorthEast');
    xlabel('Tiempo (s)');
else
    legend(['cos(\theta) Real';'C_{11} Estimada';'C_{22} Estimada'],'location','NorthEast');
    xlabel('Tiempo (s)');
    end
    subplot 312;
    hold on;
	grid;
    plot(0:Ts:300,-sin(Theta(:,2)*pi/180),'LineWidth',3);
	plot(0:Ts:300,x(6,:),'--','LineWidth',2,'Color',colors.carrot);
    xlim([0 300]);
    if(EsMatlab == 1)
    legend('-sen(\theta) Real','C_{12} Estimada','location','NorthEast');
    xlabel('Tiempo (s)');
else
    legend(['-sen(\theta) Real';'C_{12} Estimada'],'location','NorthEast');
    xlabel('Tiempo (s)');
    end

    subplot 313;
    hold on;
	grid;
    plot(0:Ts:300,sin(Theta(:,2)*pi/180),'LineWidth',3);
	plot(0:Ts:300,x(7,:),'--','LineWidth',2,'Color',colors.carrot);
    xlim([0 300]);
	suptitle('Estados c_{11},c_{12},c_{21} y c_{22}');
if(EsMatlab == 1)
    legend('sen(\theta) Real','C_{21} Estimada','location','NorthEast');
    xlabel('Tiempo (s)');
else
    legend(['sen(\theta) Real';'C_{21} Estimada'],'location','NorthEast');
    xlabel('Tiempo (s)');
    end
    
    setGraphSize(h4,wide_1);
    
if bool_print
    print('../Informe/Figuras/graf_ej5_theta','-dpdf','-bestfit');
end

	% GrÃ¡fico de correlaciÃ³n de innovaciones (debe ser ruido blanco)
	covx_g = xcorr(g(1,:)');
	covy_g = xcorr(g(2,:)');
    covvx_g = xcorr(g(3,:)');
	covvy_g = xcorr(g(4,:)');
    
	h5=figure;
    subplot 221;
	plot(covx_g)
	grid
	title('Covarianza innovaciones x')
    xlim([0 length(covx_g)]);
	
% 	figure
    subplot 222
	plot(covy_g)
	grid
	title('Covarianza innovaciones y')
    xlim([0 length(covy_g)]);
    
%     figure;
    subplot 223    
	plot(covvx_g)
	grid
	title('Covarianza innovaciones vel x')
    xlim([0 length(covvx_g)]);
	
% 	figure
    subplot 224
	plot(covvy_g)
	grid
	title('Covarianza innovaciones vel y')
	xlim([0 length(covvy_g)]);
    
    h5.Position=[0 0 1200 700];
    h5.PaperUnits='points';
    h5.PaperSize=[1200 700];
    
if bool_print
    print('../Informe/Figuras/graf_ej5_covinn','-dpdf','-bestfit');
end
	
	% Observabilidad
	Obs = obsv(Ad,C);
	rango_obs = rank(Obs);
	estados_no_observables = cant_estados - rango_obs
	
	

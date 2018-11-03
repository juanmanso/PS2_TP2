config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ7 KALMAN - Kalman Extendido
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imprimir imágenes?
bool_print=0;

acel_str = load('Acel.mat');
gyro_str = load('Gyro.mat');
radar_str = load('Radar.mat');
radar_b_v_str = load('Radar-sesgo-vel.mat');
radar_b_p_str = load('Radar-sesgo-pos.mat');
tray_str = load('trayectoria.mat');

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
cant_estados=dim*4 + 2;
cant_muestras=length(Acel);
cant_muestras_rad=length(Pradar);

% Varianzas
% var_xip=2*(0.01^3)/6;  % CHEQUEAR QUÉ VALORES PONER ACÁ
% var_xiv=0.01^2;
% var_xic=0.01^4;
% var_xib = 0.01^5;

var_xip=2*(0.01^3)/6;  % CHEQUEAR QUÉ VALORES PONER ACÁ
% var_xip=0.01^5;
var_xiv=0.01^2;
var_xic=0.01^4;
var_xib=0.01^5;

sigma_etap=10;
sigma_etav=0.1;

% Definiciones de matrices y vectores de Kalman
% Matriz de rotación
Cbe = [cos(Theta(2:2:end,2)) -sin(Theta(2:2:end,2)); sin(Theta(2:2:end,2)) cos(Theta(2:2:end,2))];

I=eye(dim); 
O=zeros(dim);

Bk1=eye(cant_estados);
C=[I O O O O;
   O I O O O];
Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,2*dim)*var_xic, var_xib, var_xib]);
R = diag([ones(1,dim)*sigma_etap^2 ones(1,dim)*sigma_etav^2]);

yk=[Pradar(:,2:3) Vradar(:,2:3)];
yk=kron(yk,[zeros(100-1,1);1]);

%Agrego sesgo a los acelerometros.
sesgo_acel_x = 0.2;
sesgo_acel_y = -0.2;

for i = 1 : cant_muestras - 1
    Acel(i,2) = Acel(i,2) + sesgo_acel_x;
    Acel(i,3) = Acel(i,3) + sesgo_acel_y;
end

% Inicio del Kalman
cov_p = [1 1]*100;
cov_v = [1 1]*0.2;
cov_c = [1 1]*1;
cov_b = [1 1]*1;

% Condiciones iniciales
x0 = [100 100 0.2 0.2 cos(40*pi/180) -sin(40*pi/180) sin(40*pi/180) cos(40*pi/180), 0, 0]';
P0_0 = diag([cov_p, cov_v, cov_c, cov_c, cov_b]);

x = x0;
P = P0_0;
xk1_k1 = x;
Pk1_k1 = P;
g = yk(1,:)';


for i=1:cant_muestras-1

    T = Ts;
    ax = Acel(i,2);
    ay = Acel(i,3);
    w = Gyro(i,2);
    bx = xk1_k1(9);
    by = xk1_k1(10);
    c11 = xk1_k1(5);
    c12 = xk1_k1(6);
    c21 = xk1_k1(7);
    c22 = xk1_k1(8);
    
    Ad = [ 1, 0, T, 0,                 (T^2*(ax - bx))/2,                 (T^2*(ay - by))/2,                                 0,                                 0,          -(T^2*c11)/2,          -(T^2*c12)/2;
           0, 1, 0, T,                                 0,                                 0,                 (T^2*(ax - bx))/2,                 (T^2*(ay - by))/2,          -(T^2*c21)/2,          -(T^2*c22)/2;
           0, 0, 1, 0, T*(ax - bx) - (T^2*w*(ay - by))/2, T*(ay - by) - (T^2*w*(ax - bx))/2,                                 0,                                 0, (c12*w*T^2)/2 - c11*T, (c11*w*T^2)/2 - c12*T;
           0, 0, 0, 1,                                 0,                                 0, T*(ax - bx) - (T^2*w*(ay - by))/2, T*(ay - by) - (T^2*w*(ax - bx))/2, (c22*w*T^2)/2 - c21*T, (c21*w*T^2)/2 - c22*T;
           0, 0, 0, 0,                   1 - (T^2*w^2)/2,                               T*w,                                 0,                                 0,                     0,                     0;
           0, 0, 0, 0,                              -T*w,                   1 - (T^2*w^2)/2,                                 0,                                 0,                     0,                     0;
           0, 0, 0, 0,                                 0,                                 0,                   1 - (T^2*w^2)/2,                               T*w,                     0,                     0;
           0, 0, 0, 0,                                 0,                                 0,                              -T*w,                   1 - (T^2*w^2)/2,                     0,                     0;
           0, 0, 0, 0,                                 0,                                 0,                                 0,                                 0,                     1,                     0;
           0, 0, 0, 0,                                 0,                                 0,                                 0,                                 0,                     0,                     1];
    
    % Predicción
    xk_k1 = Ad * xk1_k1;
    Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';

	% Corrección
    if(yk(i,1)~=0)
        gk = [innovaciones(yk(i,:),C,xk_k1)];
		Kk = Pk_k1 * C'*(R+ C*Pk_k1*C')^-1;
		xk_k = xk_k1 + Kk*(gk);
		Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
        g = [g gk];
        
        % Actualización
        xk1_k1 = xk_k;
        Pk1_k1 = Pk_k;
    else
        % Actualización
		xk1_k1 = xk_k1;
		Pk1_k1 = Pk_k1;
    end
    
%     g = [g gk];
    x = [x xk1_k1];
end

%Grafico de medida, estimada, ruidosa
h1=figure;
hold on;
grid;
plot(Preal(:,2),Preal(:,3),'LineWidth',2,'Color',[100,160,255]/255);
plot(x(1,:),x(2,:),'--','LineWidth',2,'Color',[96,24,4]/255);
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

    h1.Position=[0 0 1200 700];
    h1.PaperUnits='points';
    h1.PaperSize=[1200 700];

if bool_print
    print('../Informe/Figuras/graf_ej7','-dpdf','-bestfit');
end

	% Grafico del estado posición en función del tiempo
	h2=figure;
    subplot(3,1,[1 2]);
    hold on;
	grid;
    plot(0:Ts:300,Preal(:,2),'LineWidth',2,'Color',[100,160,255]/255);
	plot(0:Ts:300,x(1,:),'--','LineWidth',2,'Color',[96,24,4]/255);
    plot(0:Ts:300,Preal(:,3),'LineWidth',2,'Color',[249,34,74]/255);
	plot(0:Ts:300,x(2,:),'--','LineWidth',2,'Color',[16,234,38]/255);
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
    h2.Position=[0 0 700 1200];
    h2.PaperUnits='points';
    h2.PaperSize=[700 1200];
    
if bool_print
    print('../Informe/Figuras/graf_ej7_pos','-dpdf','-bestfit');
end
	
% Grafico del estado velocidad en función del tiempo
	h3=figure;
	hold on;
	grid;
    plot(0:Ts:300,Vreal(:,2),'LineWidth',2,'Color',[100,160,255]/255);
	plot(0:Ts:300,x(3,:),'--','LineWidth',2,'Color',[96,24,4]/255);
    plot(0:Ts:300,Vreal(:,3),'LineWidth',2,'Color',[249,34,74]/255);
	plot(0:Ts:300,x(4,:),'--','LineWidth',2,'Color',[16,234,38]/255);
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
    
    h3.Position=[0 0 1200 700];
    h3.PaperUnits='points';
    h3.PaperSize=[1200 700];

if bool_print
    print('../Informe/Figuras/graf_ej7_vel','-dpdf','-bestfit');
end
    
% Grafico de Cnn en función del tiempo
	h4=figure;
    subplot 311;
	hold on;
	grid;
    plot(0:Ts:300,cos(Theta(:,2)*pi/180),'LineWidth',3,'Color',[100,160,255]/255);
	plot(0:Ts:300,x(5,:),'g','LineWidth',2);
    plot(0:Ts:300,x(8,:),'--','LineWidth',2,'Color',[211,0,14]/255);
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
	plot(0:Ts:300,x(6,:),'--r','LineWidth',2);
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
	plot(0:Ts:300,x(7,:),'--r','LineWidth',2);
    xlim([0 300]);
	suptitle('Estados c_{11},c_{12},c_{21} y c_{22}');
if(EsMatlab == 1)
    legend('sen(\theta) Real','C_{21} Estimada','location','NorthEast');
    xlabel('Tiempo (s)');
else
    legend(['sen(\theta) Real';'C_{21} Estimada'],'location','NorthEast');
    xlabel('Tiempo (s)');
    end
    
    h4.Position=[0 0 1200 700];
    h4.PaperUnits='points';
    h4.PaperSize=[1200 700];
    
if bool_print
    print('../Informe/Figuras/graf_ej7_theta','-dpdf','-bestfit');
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
    print('../Informe/Figuras/graf_ej7_covinn','-dpdf','-bestfit');
end

h6=figure;
hold on;
subplot 211
plot(0:Ts:300-Ts,x(9,1:300*Fs_gyro_acel));
title('Sesgo aceleración en x');
xlabel('Tiempo (s)');
subplot 212
plot(0:Ts:300-Ts,x(10,1:300*Fs_gyro_acel));
title('Sesgo aceleración en y');
xlabel('Tiempo (s)');
	
    h6.Position=[0 0 1200 700];
    h6.PaperUnits='points';
    h6.PaperSize=[1200 700];
    
if bool_print
    print('../Informe/Figuras/graf_ej7_sesgo','-dpdf','-bestfit');
end

	% Observabilidad
	Obs = obsv(Ad,C);
	rango_obs = rank(Obs);
	estados_no_observables = cant_estados - rango_obs
	
    

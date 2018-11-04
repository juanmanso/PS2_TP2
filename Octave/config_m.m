%
% Archivo de configuraciÃ³n del main para el graph.m
%

clear all, close all;

%EsMatlab = 1;
EsMatlab = (exist('MATLAB')~=0); % Chequeo automático de MATLAB, ahora octave friendly
if(EsMatlab == 0)
    graphics_toolkit('gnuplot');
end

format long;

% Agrego las funcinoes definidas en la carpeta "funciones"
addpath('./Funciones')
addpath('./Material')

%% Definiciones para gráficos

% Colores
colors=load('colors.mat');

% Relación de aspecto
wide_1=[1200 700];
vertical_1=[700 1200];


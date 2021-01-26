%Cálculo del gradiente de presión producido por perdidas internas en la
%tubería

close all
clear all
clc

%Declaración de constantes
k_steel=0.015; %Rugosidades de los materiales
k_cu=0.001;
k_PEX=0.007;
RE=2.3438e+05;
d_int=0.0508;
dens=971.8;
v_agua=1.505;
L=100;

%Iteración para conseguir el valor de la fricción. Converge en un paso.
f_steel=zeros(1,10);
f_steel(1,1) = 0.14;

f_cu=zeros(1,10);
f_cu(1,1) = 0.14;

f_PEX=zeros(1,10);
f_PEX(1,1) = 0.14;

for i=2:10
    f_steel(1,i) = (1/(-2*log((2.51/(RE*sqrt(f_steel(1,i-1))))+(k_steel/d_int)*0.269)))^2;
    f_cu(1,i) = (1/(-2*log((2.51/(RE*sqrt(f_cu(1,i-1))))+(k_cu/d_int)*0.269)))^2;
    f_PEX(1,i) = (1/(-2*log((2.51/(RE*sqrt(f_PEX(1,i-1))))+(k_PEX/d_int)*0.269)))^2;
end

%Cálculo del pressure drop
Pd_steel = f_steel(1,10)*(L/d_int)*(dens/2)*(v_agua^2) %Pascales
Pd_cu = f_cu(1,10)*(L/d_int)*(dens/2)*(v_agua^2)
Pd_PEX = f_PEX(1,10)*(L/d_int)*(dens/2)*(v_agua^2) 
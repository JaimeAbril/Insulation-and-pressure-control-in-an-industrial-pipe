%Cálculo de los aislamientos optimos 10 años vista por minimización de
%costes

clear all
close all
clc

step=0.00001;
max=1;
e=0:step:max;
L=100;

%NOTA: LOS CÁLCULOS DE LOS COSTES DE LA ENERGÍA ESTAN POR METRO DE
%CONDUCCIÓN

%Constantes de los materiales

k_steel = 50; %W/m*k
k_cu = 364;
k_PEX = 0.36;

%Listado de precios para los aislantes por milímetro de espesor y por metro
%longitudinal

P_fv = 0.147;
P_foam = 0.828;
P_PE = 0.064;

%Listado de precios de los materiales para la tubería por metro
%longitudinal

P_steel = 7.00;
P_cu = 18.79;
P_PEX = 2.20;

%Coste fijo para la instalación (por si queremos fijar algun valor de mano de obra o similares)

Cf = 0;

%% Cálculos previos

%Datos a tempreatura media de operación 80ºC

vel_agua = 1.505; %(m/s)
vel_aire = 5.56;
d_agua = 971.8; %(kg/m3)
d_aire = 1;
visd_agua = 0.317*10^-3; %(N*s/m^2)
visd_aire = 2.09*10^-5;
Pr_agua = 2.22; %Número de Prandt
k_agua = 0.673; %Conductividad del agua (W/m*ºC)
D_int = 0.0508; %Diametro interior de la tubería
D_ext = 0.0603; %Diametro exterior de la tubería

RE = (d_agua*vel_agua*D_int)/visd_agua; %Número de Reynolds
Nu = 0.0269 * RE^0.8 * Pr_agua^(1/3); %Número de Nusselt
h_agua = (Nu * k_agua)/D_int; %Coefficiente de transmisión de calor del agua por convección

RE_aire = (d_aire*vel_aire*D_ext)/visd_aire; %Número de Reynolds para el aire

%% Acero

Rt_steel = (1/(4*pi*(D_int/2)^2*h_agua))+((log((D_ext/2)/(D_int/2)))/(2*pi*L*k_steel))+(1/(293.58+9424.7*0));
Rt_steel_fv = zeros(1,length(e));
Rt_steel_foam = zeros(1,length(e));
Rt_steel_PE = zeros(1,length(e));

Q_steel = (90-20)/Rt_steel;
Q_steel_fv = zeros(1,length(e));
Q_steel_foam = zeros(1,length(e));
Q_steel_PE = zeros(1,length(e));

Ce_steel = (Q_steel/1000)*0,19*10*365*24*100;
Ce_steel_fv = zeros(1,length(e));
Ce_steel_foam = zeros(1,length(e));
Ce_steel_PE = zeros(1,length(e));

Ci_steel = Cf+P_fv*0*100*1000+100*P_steel;
Ci_steel_fv = zeros(1,length(e));
Ci_steel_foam = zeros(1,length(e));
Ci_steel_PE = zeros(1,length(e));

Ct_steel = Ce_steel + Ci_steel;
Ct_steel_fv = zeros(1,length(e));
Ct_steel_foam = zeros(1,length(e));
Ct_steel_PE = zeros(1,length(e));

for i=1:length(e)
    Rt_steel_fv(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_steel))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.035))+(log((e(1,i)+(0.03115))/(e(1,i)+0.03015))/245735.37)+(1/(293.58+9424.7*e(1,i)));
    Rt_steel_foam(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_steel))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.042))+(1/(293.58+9424.7*e(1,i)));
    Rt_steel_PE(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_steel))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.035))+(1/(293.58+9424.7*e(1,i)));
end

for i=1:length(e)
    Q_steel_fv(1,i) = ((90-20)/(Rt_steel_fv(1,i)));
    Q_steel_foam(1,i) = (90-20)/Rt_steel_foam(1,i);
    Q_steel_PE(1,i) = (90-20)/Rt_steel_PE(1,i);
end

for i=1:length(e)
    Ce_steel_fv(1,i) = (Q_steel_fv(1,i)/1000)*0.19*10*365*24;
    Ce_steel_foam(1,i) = (Q_steel_foam(1,i)/1000)*0.19*10*365*24;
    Ce_steel_PE(1,i) = (Q_steel_PE(1,i)/1000)*0.19*10*365*24;
end

for i=1:length(e)
    Ci_steel_fv(1,i) = Cf+P_fv*e(1,i)*100*1000+100*P_steel;
    Ci_steel_foam(1,i) = Cf+P_foam*e(1,i)*100*1000+100*P_steel;
    Ci_steel_PE(1,i) = Cf+P_PE*e(1,i)*100*1000+100*P_steel;
end

for i=1:length(e)
    Ct_steel_fv(1,i) = Ce_steel_fv(1,i)+Ci_steel_fv(1,i);
    Ct_steel_foam(1,i) = Ce_steel_foam(1,i)+Ci_steel_foam(1,i);
    Ct_steel_PE(1,i) = Ce_steel_PE(1,i)+Ci_steel_PE(1,i);
end


%% Cobre

Rt_cu = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_cu))+(1/(293.58+9424.7*0));
Rt_cu_fv = zeros(1,length(e));
Rt_cu_foam = zeros(1,length(e));
Rt_cu_PE = zeros(1,length(e));

Q_cu = (90-20)/Rt_cu;
Q_cu_fv = zeros(1,length(e));
Q_cu_foam = zeros(1,length(e));
Q_cu_PE = zeros(1,length(e));

Ce_cu = (Q_cu/1000)*0,19*10*365*24*100;
Ce_cu_fv = zeros(1,length(e));
Ce_cu_foam = zeros(1,length(e));
Ce_cu_PE = zeros(1,length(e));

Ci_cu = Cf+P_fv*0*100*1000+100*P_cu;
Ci_cu_fv = zeros(1,length(e));
Ci_cu_foam = zeros(1,length(e));
Ci_cu_PE = zeros(1,length(e));

Ct_cu = Ce_cu + Ci_cu;
Ct_cu_fv = zeros(1,length(e));
Ct_cu_foam = zeros(1,length(e));
Ct_cu_PE = zeros(1,length(e));

for i=1:length(e)
    Rt_cu_fv(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_cu))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.035))+(log((e(1,i)+(0.03115))/(e(1,i)+0.03015))/245735.37)+(1/(293.58+9424.7*e(1,i)));
    Rt_cu_foam(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_cu))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.042))+(1/(293.58+9424.7*e(1,i)));
    Rt_cu_PE(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_cu))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.035))+(1/(293.58+9424.7*e(1,i)));
end

for i=1:length(e)
    Q_cu_fv(1,i) = (90-20)/Rt_cu_fv(1,i);
    Q_cu_foam(1,i) = (90-20)/Rt_cu_foam(1,i);
    Q_cu_PE(1,i) = (90-20)/Rt_cu_PE(1,i);
end

for i=1:length(e)
    Ce_cu_fv(1,i) = (Q_cu_fv(1,i)/1000)*0.19*10*365*24;
    Ce_cu_foam(1,i) = (Q_cu_foam(1,i)/1000)*0.19*10*365*24;
    Ce_cu_PE(1,i) = (Q_cu_PE(1,i)/1000)*0.19*10*365*24;
end

for i=1:length(e)
    Ci_cu_fv(1,i) = Cf+P_fv*e(1,i)*100*1000+100*P_cu;
    Ci_cu_foam(1,i) = Cf+P_foam*e(1,i)*100*1000+100*P_cu;
    Ci_cu_PE(1,i) = Cf+P_PE*e(1,i)*100*1000+100*P_cu;
end

for i=1:length(e)
    Ct_cu_fv(1,i) = Ce_cu_fv(1,i)+Ci_cu_fv(1,i);
    Ct_cu_foam(1,i) = Ce_cu_foam(1,i)+Ci_cu_foam(1,i);
    Ct_cu_PE(1,i) = Ce_cu_PE(1,i)+Ci_cu_PE(1,i);
end


%% PEX

Rt_PEX = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_PEX))+(1/(293.58+9424.7*0));
Rt_PEX_fv = zeros(1,length(e));
Rt_PEX_foam = zeros(1,length(e));
Rt_PEX_PE = zeros(1,length(e));

Q_PEX = (90-20)/Rt_PEX;
Q_PEX_fv = zeros(1,length(e));
Q_PEX_foam = zeros(1,length(e));
Q_PEX_PE = zeros(1,length(e));

Ce_PEX = (Q_PEX/1000)*0,19*10*365*24*100;
Ce_PEX_fv = zeros(1,length(e));
Ce_PEX_foam = zeros(1,length(e));
Ce_PEX_PE = zeros(1,length(e));

Ci_PEX = Cf+P_fv*0*100*1000+100*P_PEX;
Ci_PEX_fv = zeros(1,length(e));
Ci_PEX_foam = zeros(1,length(e));
Ci_PEX_PE = zeros(1,length(e));

Ct_PEX = Ce_PEX + Ci_PEX;
Ct_PEX_fv = zeros(1,length(e));
Ct_PEX_foam = zeros(1,length(e));
Ct_PEX_PE = zeros(1,length(e));

for i=1:length(e)
    Rt_PEX_fv(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_PEX))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.035))+(log((e(1,i)+(0.03115))/(e(1,i)+0.03015))/245735.37)+(1/(293.58+9424.7*e(1,i)));
    Rt_PEX_foam(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_PEX))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.042))+(1/(293.58+9424.7*e(1,i)));
    Rt_PEX_PE(1,i) = (1/(4*pi*(D_int/2)^2*h_agua))+(log((D_ext/2)/(D_int/2))/(2*pi*L*k_PEX))+((log(1+(e(1,i)/0.03015)))/(2*pi*100*0.035))+(1/(293.58+9424.7*e(1,i)));
end

for i=1:length(e)
    Q_PEX_fv(1,i) = (90-20)/Rt_PEX_fv(1,i);
    Q_PEX_foam(1,i) = (90-20)/Rt_PEX_foam(1,i);
    Q_PEX_PE(1,i) = (90-20)/Rt_PEX_PE(1,i);
end

for i=1:length(e)
    Ce_PEX_fv(1,i) = (Q_PEX_fv(1,i)/1000)*0.19*10*365*24;
    Ce_PEX_foam(1,i) = (Q_PEX_foam(1,i)/1000)*0.19*10*365*24;
    Ce_PEX_PE(1,i) = (Q_PEX_PE(1,i)/1000)*0.19*10*365*24;
end

for i=1:length(e)
    Ci_PEX_fv(1,i) = Cf+P_fv*e(1,i)*100*1000+100*P_PEX;
    Ci_PEX_foam(1,i) = Cf+P_foam*e(1,i)*100*1000+100*P_PEX;
    Ci_PEX_PE(1,i) = Cf+P_PE*e(1,i)*100*1000+100*P_PEX;
end

for i=1:length(e)
    Ct_PEX_fv(1,i) = Ce_PEX_fv(1,i)+Ci_PEX_fv(1,i);
    Ct_PEX_foam(1,i) = Ce_PEX_foam(1,i)+Ci_PEX_foam(1,i);
    Ct_PEX_PE(1,i) = Ce_PEX_PE(1,i)+Ci_PEX_PE(1,i);
end


%% Plotting

%Plot para el acero
figure
plot(e,Ct_steel_fv,'LineWidth',1.5)
title('Total cost for the installation made of steel', 'Fontsize', 21)
set(gca,'FontSize',20)
xlabel('Thickness (m)', 'Fontsize', 25)
ylabel('Cost per meter of conduction along 10 years operating', 'Fontsize', 15)
grid on
hold on
plot(e,Ct_steel_foam,'r','LineWidth',1.5)
hold on
plot(e,Ct_steel_PE,'LineWidth',1)
% xline(0.11, '--')
% yline(0.1663,'--')
% yline(3.4801,'--')
% yline(1.4438,'--')
legend({'Fiberglass insulation','Elastomeric foam', 'Polyethylene'},'Location','northeast', 'Fontsize', 20)
hold off

%Plot para el cobre
figure
plot(e,Ct_cu_fv,'LineWidth',1.5)
title('Total cost for the installation made of cooper', 'Fontsize', 21)
set(gca,'FontSize',20)
xlabel('Thickness (m)', 'Fontsize', 25)
ylabel('Cost per meter of conduction along 10 years operating', 'Fontsize', 15)
grid on
hold on
plot(e,Ct_cu_foam,'r','LineWidth',1.5)
hold on
plot(e,Ct_cu_PE,'LineWidth',1)
legend({'Fiberglass insulation','Elastomeric foam', 'Polyethylene'},'Location','northeast', 'Fontsize', 20)
hold off


%Plot para el PEX
figure
plot(e,Ct_PEX_fv,'LineWidth',1.5)
title('Total cost for the installation made of Cross-linked Polyethylene', 'Fontsize', 21)
set(gca,'FontSize',20)
xlabel('Thickness (m)', 'Fontsize', 25)
ylabel('Cost per meter of conduction along 10 years operating', 'Fontsize', 15)
grid on
hold on
plot(e,Ct_steel_foam,'r','LineWidth',1.5)
hold on
plot(e,Ct_steel_PE,'LineWidth',1)
legend({'Fiberglass insulation','Elastomeric foam', 'Polyethylene'},'Location','northeast', 'Fontsize', 20)
hold off

%Plot solo PE -- NO GRAPH -- SE DERRITE
figure
plot(e,Ct_steel_PE,'LineWidth',1.5)
title('Total cost for the installation insulated in PE', 'Fontsize', 21)
set(gca,'FontSize',20)
xlabel('Thickness (m)', 'Fontsize', 25)
ylabel('Cost per meter of conduction along 10 years operating', 'Fontsize', 15)
grid on
hold on
plot(e,Ct_cu_PE,'r','LineWidth',1.5)
hold on
plot(e,Ct_PEX_PE,'LineWidth',1)
legend({'Steel','Cu', 'PEX'},'Location','northeast', 'Fontsize', 20)
hold off

%Plot para la espuma
figure
plot(e,Ct_steel_foam,'LineWidth',1.5)
title('Total cost for the installation insulated in Elastomeric Foam', 'Fontsize', 21)
set(gca,'FontSize',20)
xlabel('Thickness (m)', 'Fontsize', 25)
ylabel('Cost per meter of conduction along 10 years operating', 'Fontsize', 15)
grid on
hold on
plot(e,Ct_cu_foam,'r','LineWidth',1.5)
hold on
plot(e,Ct_PEX_foam,'LineWidth',1)
legend({'Steel','Cu', 'PEX'},'Location','northeast', 'Fontsize', 20)
hold off

%Plot para la fibra de vidrio
figure
plot(e,Ct_steel_fv,'LineWidth',1.5)
title('Total cost for the installation insulated in Fiberglass', 'Fontsize', 21)
set(gca,'FontSize',20)
xlabel('Thickness (m)', 'Fontsize', 25)
ylabel('Cost per meter of conduction along 10 years operating', 'Fontsize', 15)
grid on
hold on
plot(e,Ct_cu_fv,'r','LineWidth',1.5)
hold on
plot(e,Ct_PEX_fv,'LineWidth',1)
legend({'Steel','Cu', 'PEX'},'Location','northeast', 'Fontsize', 20)
hold off

%Comparacion steel/cu-foam/fv
figure
plot(e,Ct_steel_fv,'LineWidth',1.5)
title('Total cost comparison for Steel and Cu', 'Fontsize', 21)
set(gca,'FontSize',20)
xlabel('Thickness (m)', 'Fontsize', 25)
ylabel('Cost per meter of conduction along 10 years operating', 'Fontsize', 15)
grid on
hold on
plot(e,Ct_cu_fv,'r','LineWidth',1.5)
hold on
plot(e,Ct_steel_foam,'LineWidth',1)
hold on
plot(e,Ct_cu_foam,'LineWidth',1.5)
xline(0.2440,'--')
xline(0.0885,'--')
legend({'Steel & Fiberglass','Cu & Fiberglass', 'Steel & Elastomeric Foam', 'Cu & Elastomeric Foam'},'Location','northeast', 'Fontsize', 20)
hold off

%% Cálculo mínimos

[a1,b1]=min(Ct_steel_foam);
e_opt_steel_foam=e(1,b1); %--->0.0885

[a2,b2]=min(Ct_cu_foam);
e_opt_cu_foam=e(1,b2); %--->0.0886

[a3,b3]=min(Ct_steel_fv);
e_opt_steel_fv=e(1,b3); %--->0.2439

[a4,b4]=min(Ct_cu_fv);
e_opt_cu_fv=e(1,b4); %--->0.2440

incremento_precio_cu_foam=(a2-a1)/a1*100; %---> 4.6069 por ciento respecto al steel + ganamos menos rugosidad (motor barato)
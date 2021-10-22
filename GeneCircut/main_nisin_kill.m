%{
global T_0_1 gamma_M_1 alpha_c  gamma_c K_1 alpha_a  gamma_a %K_a K_b K_ia 
global gamma_M_2 gamma_M_3
global alpha_D  gamma_D alpha_N gamma_N alpha_L gamma_L K_AP K_A m mu K_two Kapi kapi K_b

T_0_1=0.1; %basal transcription rate
gamma_M_1=0;
alpha_c=1; %effective factor of AgrC synthesis
gamma_c=2;%degradation rate of AgrC
K_1=1; %assosciation rate of AgrC&AIP1
alpha_a=1;%effective factor of AgrA synthesis
gamma_a=2; %degradation rate of unphosphorylated AgrA
gamma_M_2=2;%degradation rate of RNA2
gamma_M_3=2;%degradation rate of RNA3
alpha_D=1;
gamma_D=2;%degradation rate of DNAsei
alpha_N=1;
gamma_N=2;%degradation rate of NisinPV
alpha_L=1;
gamma_L=2;%degradation rate of lysisE7
K_AP=1;%rate of dephosphorylation of AgrA
K_A=10;%rate of phosphorylation of AgrA
m=50;
Kapi=10; %Maximal AgrA dependent transcription rate
kapi=1; %conc. 0f AgrA_P for half maximal transcription rate
K_b=0.1; %basal transcription rate of P2
mu=1; %volume of Ecoli u/m^3
K_two=0.1;%dissosciation rate of AgrC&AIP1;

tlast = 20 ; 
figure
hold on

AgrC = 0.1;
AIP_1=0.1;
AgrA = 0.1;
AgrC_AIP = 0;
AgrA_P = 0;
M_2 = 0;
M_3 = 0;
DNAse_i = 0;
Nisin_PV = 0;
Lysis_E7 = 0;

statevar_i = [AgrC; AIP_1; AgrA; AgrC_AIP; AgrA_P; M_2; M_3; DNAse_i; Nisin_PV; Lysis_E7];
[time,statevars] = ode45(@dydt_kill,[0,tlast],statevar_i) ;

AgrC = statevars(:,1);
AIP_1 = statevars(:,2);
AgrA = statevars(:,3);
AgrC_AIP = statevars(:,4);
AgrA_P = statevars(:,5);
M_2 = statevars(:,6);
M_3 = statevars(:,7);
DNAse_i = statevars(:,8);
Nisin_PV = statevars(:,9);
Lysis_E7 = statevars(:,10);

plot(time,Nisin_PV,'r','LineWidth',2.25)
xlabel('time')
ylabel('Internal NisinPV concentration')

figure
hold on
plot(time,Lysis_E7,'g','LineWidth',2.25)
xlabel('time')
ylabel('Internal LysisE7 concentration')
%}
close all
clear all
clc

global alpha_c m Kapi v kapi K_b gamma_c K_1 K_two alpha_a gamma_a K_AP K_A
global gamma_M_2 gamma_M_3 alpha_D gamma_D alpha_N gamma_N alpha_L gamma_L D_1 D_2

alpha_c = 1; %effective factor of AgrC synthesis
m = 1;
Kapi = 10; %Maximal AgrA dependent transcription rate
v = 1;
kapi = 1; %conc. 0f AgrA_P for half maximal transcription rate
K_b = 0.1; %basal transcription rate of P2
gamma_c = 2; %degradation rate of AgrC
K_1 = 1; %assosciation rate of AgrC&AIP1
K_two = 0.1; %dissosciation rate of AgrC&AIP1;
alpha_a = 1; %effective factor of AgrA synthesis
gamma_a = 2; %degradation rate of unphosphorylated AgrA
K_AP = 1; %rate of dephosphorylation of AgrA
K_A = 10; %rate of phosphorylation of AgrA
gamma_M_2 = 2; %degradation rate of RNA2
gamma_M_3 = 2; %degradation rate of RNA3
alpha_D = 1;
gamma_D = 2; %degradation rate of DNAsei
alpha_N = 1;
gamma_N = 2; %degradation rate of NisinPV
alpha_L = 0.1;
gamma_L = 4; %degradation rate of lysisE7
D_1 = 2;%degradation rateof AgrC&AIP1
D_2 = 2;%degradation rateof AgrA_P 

tlast = 20 ; 

AIP_arr=0.01:0.001:10;

AgrC = 0.1;
AIP_1=0.01;
AgrA = 0.1;
AgrC_AIP = 0;
AgrA_P = 0;
M_2 = 0;
M_3 = 0;
DNAse_i = 0;
Nisin_PV = 0;
Lysis_E7 = 0;

statevar_i = [AgrC; AIP_1; AgrA; AgrC_AIP; AgrA_P; M_2; M_3; DNAse_i; Nisin_PV; Lysis_E7];
[time,statevars] = ode45(@dydt_kill,[0,tlast],statevar_i) ;

AgrC = statevars(:,1);
AIP_1 = statevars(:,2);
AgrA = statevars(:,3);
AgrC_AIP = statevars(:,4);
AgrA_P = statevars(:,5);
M_2 = statevars(:,6);
M_3 = statevars(:,7);
DNAse_i = statevars(:,8);
Nisin_PV = statevars(:,9);
Lysis_E7 = statevars(:,10);
figure
plot(time,Nisin_PV,'r','LineWidth',2.5)
title('NisinPV concentration v/s time','FontSize',20)
xlabel('time(hr)','FontSize',20)
ylabel(' NisinPV concentration(\mu mol)','FontSize',20)
figure
plot(AIP_1,Nisin_PV,'y','LineWidth',2.5)
title('NisinPV conc v/s AIP-1 conc','FontSize',20)
xlabel(' AIP1 concentration','FontSize',20)
ylabel(' NisinPV concentration(\mu mol)','FontSize',20)
figure
hold on
plot(time,Lysis_E7,'g','LineWidth',2.5)
title('Lysis E7 concentration v/s time','FontSize',20)
xlabel('time(hr)','FontSize',20)
ylabel(' concentration of LysisE7(\mu mol)','FontSize',20)
hold off
figure
hold on
plot(time,AIP_1,'m','LineWidth',2.5)
title('AIP-1 concentration v/s time','FontSize',20)
xlabel('time(hr)','FontSize',20)
ylabel(' concentration of AIP-1(\mu mol)','FontSize',20)

figure
hold on
plot(time,DNAse_i,'b','LineWidth',2.5)
title('DNAse-i concentration v/s time','FontSize',20)
xlabel('time(hr)','FontSize',20)
ylabel('DNAse-i concentration(\mu mol)','FontSize',20)
figure
plot(time,log(Nisin_PV),'r','LineWidth',2.5)
title('log of NisinPV conc. v/s time','FontSize',20)
xlabel('time(hr)','FontSize',20)
ylabel('log of NisinPV concentration','FontSize',20)
figure
hold on
plot(time,log(Lysis_E7),'g','LineWidth',2.5)
plot(time,log(DNAse_i),'b','LineWidth',2.5)
title('log of conc. v/s time','FontSize',20)
legend('log(LysisE7)','log(DNAsei)')
xlabel('time(hr)','FontSize',20)
ylabel('log of concentration' ,'FontSize',20)
hold off
figure
hold on
plot(time,log(AIP_1),'m','LineWidth',2.5)
title('log of AIP-1 conc v/s time','FontSize',20)
xlabel('time(hr)','FontSize',20)
ylabel(' log of conc of AIP-1','FontSize',20)



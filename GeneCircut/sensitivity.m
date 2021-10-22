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


AIP_arr=0.01:0.01:1;
%Nisin_max=zeros(length(AIP_arr));
Nisin_max=[];

for i=1:length(AIP_arr)
    AgrC = 0.1;
    AIP_1=AIP_arr(i);
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
    %Nisin_max(i)=max(Nisin_PV);
    Nisin_max=[Nisin_max,max(Nisin_PV)];
   
    
    
end
figure
scatter(AIP_arr,Nisin_max)
title('Sensitivity plot of NisinPV','FontSize',20)
xlabel('Initial AIP conc. (\mu mol)','FontSize',20)
ylabel(' Max nisin PV conc.(\mu mol)','FontSize',20)
figure
scatter(AIP_arr,log(Nisin_max))
title('Sensitivity plot of log of NisinPV','FontSize',20)
xlabel('Initial AIP conc. (\mu mol)','FontSize',20)
ylabel(' log of Max nisin PV conc','FontSize',20)



%{
function deriv = dydt_kill(t,statevar)

global alpha_c gamma_c K_1 alpha_a gamma_a %K_a K_b K_ia
global gamma_M_2 gamma_M_3
global alpha_D  gamma_D alpha_N  gamma_N alpha_L gamma_L K_AP K_A m v K_two Kapi kapi K_b

AgrC = statevar(1);
AIP_1 = statevar(2);
AgrA = statevar(3);
AgrC_AIP = statevar(4);
AgrA_P = statevar(5);
M_2 = statevar(6);
M_3 = statevar(7);
DNAse_i = statevar(8);
Nisin_PV = statevar(9);
Lysis_E7 = statevar(10);


dAgrC = (alpha_c*m*Kapi*((AgrA_P./(m*v))./(kapi+(AgrA_P./(m*v)))))+ alpha_c*m*K_b - gamma_c*AgrC - K_1*AgrC*AIP_1 + K_two*AgrC_AIP ;

dAIP_1 = -K_1*AgrC*AIP_1 + K_two*AgrC_AIP;

dAgrA = alpha_a*m*Kapi*((AgrA_P./(m*v))./(kapi+(AgrA_P./(m*v))))+ alpha_a*m*K_b - gamma_a*AgrA + K_AP*AgrA_P - (K_A/(m*v))*AgrA*AgrC_AIP;
%{- (K_a*K_b*AgrC_AIP*AgrA*P)/(K_ia*K_b+K_b*AgrA+K_a*P+AgrA*P)} 

dAgrC_AIP = K_1*AgrC*AIP_1 - (K_A/(m*v))*AgrA*AgrC_AIP - K_two*AgrC_AIP; 
%(K_a*K_b*AgrC_AIP*AgrA*P)/(K_ia*K_b+K_b*AgrA+K_a*P+AgrA*P);

dAgrA_P = (K_A/(m*v))*AgrA*AgrC_AIP -K_AP*AgrA_P; 
%- (K_a*K_b*AgrC_AIP*AgrA*P)/(K_ia*K_b+K_b*AgrA+K_a*P+AgrA*P);

dM_2 = m*Kapi*((AgrA_P./(m*v))./(kapi+(AgrA_P./(m*v))))+m*K_b  - gamma_M_2*M_2 ; 
%- alpha_D*exp(-u_2*T_D)*M2(t-T_D) -alpha_N*exp(-u_2*T_N)*M2(t-T_N);

dM_3 = m*Kapi*((AgrA_P./(m*v))./(kapi+(AgrA_P./(m*v)))) - gamma_M_3*M_3 ; 
%- alpha_L*exp(-u_3*T_L)*M3(t-T_L);

dDNAse_i = alpha_D*m*Kapi*((AgrA_P./(m*v))./(kapi+(AgrA_P./(m*v)))) - gamma_D*DNAse_i ;

dNisin_PV = alpha_N*m*Kapi*((AgrA_P./(m*v))./(kapi+(AgrA_P./(m*v)))) - gamma_N*Nisin_PV ;

dLysis_E7 = alpha_L*m*Kapi*((AgrA_P./(m*v))./(kapi+(AgrA_P./(m*v)))) - gamma_L*Lysis_E7;


deriv = [dAgrC; dAIP_1; dAgrA; dAgrC_AIP; dAgrA_P; dM_2; dM_3; dDNAse_i; dNisin_PV; dLysis_E7];

return
%}

function deriv = dydt_kill(t,statevar)

global alpha_c m Kapi v kapi K_b gamma_c K_1 K_two alpha_a gamma_a K_AP K_A
global gamma_M_2 gamma_M_3 alpha_D gamma_D alpha_N gamma_N alpha_L gamma_L D_1 D_2

AgrC = statevar(1);
AIP_1 = statevar(2);
AgrA = statevar(3);
AgrC_AIP = statevar(4);
AgrA_P = statevar(5);
M_2 = statevar(6);
M_3 = statevar(7);
DNAse_i = statevar(8);
Nisin_PV = statevar(9);
Lysis_E7 = statevar(10);

dAgrC = alpha_c*m*Kapi*((AgrA_P/(m*v))/(kapi+(AgrA_P/(m*v))))+ alpha_c*m*K_b - gamma_c*AgrC - K_1*AgrC*AIP_1 + K_two*AgrC_AIP;

dAIP_1 = -K_1*AgrC*AIP_1 + K_two*AgrC_AIP;

dAgrA = alpha_a*m*Kapi*((AgrA_P/(m*v))/(kapi+(AgrA_P/(m*v))))+ alpha_a*m*K_b - gamma_a*AgrA + K_AP*AgrA_P - (K_A/(m*v))*AgrA*AgrC_AIP;

dAgrC_AIP = K_1*AgrC*AIP_1 - (K_A/(m*v))*AgrA*AgrC_AIP - K_two*AgrC_AIP-D_1*AgrC_AIP;

dAgrA_P = (K_A/(m*v))*AgrA*AgrC_AIP -K_AP*AgrA_P-D_2*AgrA_P;

dM_2 = m*Kapi*((AgrA_P/(m*v))/(kapi+(AgrA_P/(m*v))))+m*K_b  - gamma_M_2*M_2;

dM_3 = m*Kapi*((AgrA_P/(m*v))/(kapi+(AgrA_P/(m*v)))) - gamma_M_3*M_3;

dDNAse_i = alpha_D*m*Kapi*((AgrA_P/(m*v))/(kapi+(AgrA_P/(m*v)))) - gamma_D*DNAse_i;

dNisin_PV = alpha_N*m*Kapi*((AgrA_P/(m*v))/(kapi+(AgrA_P/(m*v)))) - gamma_N*Nisin_PV;

dLysis_E7 = alpha_L*m*Kapi*((AgrA_P/(m*v))/(kapi+(AgrA_P/(m*v)))) - gamma_L*Lysis_E7;

deriv = [dAgrC; dAIP_1; dAgrA; dAgrC_AIP; dAgrA_P; dM_2; dM_3; dDNAse_i; dNisin_PV; dLysis_E7];

return
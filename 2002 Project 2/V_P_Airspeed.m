% Andrew Logue
% House Cleaning
clear all
clc
close
load('VenturiDataFile.mat')
load('PitotProbeDataFile.mat')
%%
% Base Calculations
N = height(VV1); % Number of elements 
I = 500; 
 VV1 = table2array(VV1);
 VV2 = table2array(VV2);
 VV3 = table2array(VV3);
 VV4 = table2array(VV4);
 VV5 = table2array(VV5);
 VV6 = table2array(VV6);
 VV7 = table2array(VV7);
 VV8 = table2array(VV8);
 VV9 = table2array(VV9);
 VV10 = table2array(VV10);
 VV11 = table2array(VV11);
 VV12 = table2array(VV12);
%average value of the atmospheric temperature and pressure 
Patm = (sum(VV1(:,1)) + sum(VV2(:,1)) + sum(VV3(:,1)) + sum(VV4(:,1)) + sum(VV5(:,1)) + sum(VV6(:,1)) + sum(VV7(:,1)) + sum(VV8(:,1)) + sum(VV9(:,1)) + sum(VV10(:,1)) + sum(VV11(:,1)) + sum(VV12(:,1))) / (N*12);
Tatm = (sum(VV1(:,2)) + sum(VV2(:,2)) + sum(VV3(:,2)) + sum(VV4(:,2)) + sum(VV5(:,2)) + sum(VV6(:,2)) + sum(VV7(:,2)) + sum(VV8(:,2)) + sum(VV9(:,2)) + sum(VV10(:,2)) + sum(VV11(:,2)) + sum(VV12(:,2))) / (N*12);
R = 8.314; %[J/mol·K]

ratio = 1/9.5;
PP1 = table2array(PP1);
 PP2 = table2array(PP2);
 PP3 = table2array(PP3);
 PP4 = table2array(PP4);
 PP5 = table2array(PP5);
 PP6 = table2array(PP6);
 PP7 = table2array(PP7);
 PP8 = table2array(PP8);
 PP9 = table2array(PP9);
 PP10 = table2array(PP10);
 PP11 = table2array(PP11);
 PP12 = table2array(PP12);
%average value of the atmospheric temperature and pressure 
Patm2 = (sum(PP1(:,1)) + sum(PP2(:,1)) + sum(PP3(:,1)) + sum(PP4(:,1)) + sum(PP5(:,1)) + sum(PP6(:,1)) + sum(PP7(:,1)) + sum(PP8(:,1)) + sum(PP9(:,1)) + sum(PP9(:,1)) + sum(PP11(:,1)) + sum(PP12(:,1))) / (N*12);
Tatm2 = (sum(PP1(:,2)) + sum(PP2(:,2)) + sum(PP3(:,2)) + sum(PP4(:,2)) + sum(PP5(:,2)) + sum(PP6(:,2)) + sum(PP7(:,2)) + sum(PP8(:,2)) + sum(PP9(:,2)) + sum(PP9(:,2)) + sum(PP11(:,2)) + sum(PP12(:,2))) / (N*12);

%%
% Obtaining the various pressures with respect to their voltage
for Q = 1:500
	%VOLTAGE = 0.5
    VD5(:,1) = VV3(Q,3) + VV4(Q,3) + VV7(Q,3) + VV8(Q,3) + VV11(Q,3) + VV12(Q,3);
    %VD5(:,2) = VV3(Q,4) + VV4(Q,4) + VV7(Q,4) + VV8(Q,4) + VV11(Q,4) + VV12(Q,4);
    %VOLTAGE = 2
    V2(:,1) = VV1(Q,3) + VV2(Q,3) + VV5(Q,3) + VV6(Q,3) + VV9(Q,3) + VV10(Q,3);
    %V2(:,2) = VV1(Q,3) + VV2(Q,4) + VV5(Q,4) + VV6(Q,4) + VV9(Q,4) + VV10(Q,4);
    
    %VOLTAGE = 2.5
    V2D5(:,1) = VV3(Q+I,3) + VV4(Q+I,3) + VV7(Q+I,3) + VV8(Q+I,3) + VV11(Q+I,3) + VV12(Q+I,3);
    %V2D5(:,2) = VV3(Q+I,4) + VV4(Q+I,4) + VV7(Q+I,4) + VV8(Q+I,4) + VV11(Q+I,4) + VV12(Q+I,4);
    %VOLTAGE = 4
    V4(:,1) = VV1(Q+I,3) + VV2(Q+I,3) + VV5(Q+I,3) + VV6(Q+I,3) + VV9(Q+I,3) + VV10(Q+I,3);
    %V4(:,2) = VV1(Q+I,3) + VV2(Q+I,4) + VV5(Q+I,4) + VV6(Q+I,4) + VV9(Q+I,4) + VV10(Q+I,4);
   
    %VOLTAGE = 4.5
    V4D5(:,1) = VV3(Q+2*I,3) + VV4(Q+2*I,3) + VV7(Q+2*I,3) + VV8(Q+2*I,3) + VV11(Q+2*I,3) + VV12(Q+2*I,3);
    %V4D5(:,2) = VV3(Q+2*I,4) + VV4(Q+2*I,4) + VV7(Q+2*I,4) + VV8(Q+2*I,4) + VV11(Q+2*I,4) + VV12(Q+2*I,4);
    %VOLTAGE = 6
    V6(:,1) = VV1(Q+2*I,3) + VV2(Q+2*I,3) + VV5(Q+2*I,3) + VV6(Q+2*I,3) + VV9(Q+2*I,3) + VV10(Q+2*I,3);
    %V6(:,2) = VV1(Q+2*I,3) + VV2(Q+2*I,4) + VV5(Q+2*I,4) + VV6(Q+2*I,4) + VV9(Q+2*I,4) + VV10(Q+2*I,4);

    %VOLTAGE = 6.5
    V6D5(:,1) = VV3(Q+3*I,3) + VV4(Q+3*I,3) + VV7(Q+3*I,3) + VV8(Q+3*I,3) + VV11(Q+3*I,3) + VV12(Q+3*I,3);
    %V6D5(:,2) = VV3(Q+3*I,4) + VV4(Q+3*I,4) + VV7(Q+3*I,4) + VV8(Q+3*I,4) + VV11(Q+3*I,4) + VV12(Q+3*I,4);
    %VOLTAGE = 8
    V8(:,1) = VV1(Q+3*I,3) + VV2(Q+3*I,3) + VV5(Q+3*I,3) + VV6(Q+3*I,3) + VV9(Q+3*I,3) + VV10(Q+3*I,3);
    %V8(:,2) = VV1(Q+3*I,3) + VV2(Q+3*I,4) + VV5(Q+3*I,4) + VV6(Q+3*I,4) + VV9(Q+3*I,4) + VV10(Q+3*I,4);

    %VOLTAGE = 8.5
    V8D5(:,1) = VV3(Q+4*I,3) + VV4(Q+4*I,3) + VV7(Q+4*I,3) + VV8(Q+4*I,3) + VV11(Q+4*I,3) + VV12(Q+4*I,3);
    %V8D5(:,2) = VV3(Q+4*I,4) + VV4(Q+4*I,4) + VV7(Q+4*I,4) + VV8(Q+4*I,4) + VV11(Q+4*I,4) + VV12(Q+4*I,4);
    %VOLTAGE = 10
    V10(:,1) = VV1(Q+4*I,3) + VV2(Q+4*I,3) + VV5(Q+4*I,3) + VV6(Q+4*I,3) + VV9(Q+4*I,3) + VV10(Q+4*I,3);
    %V10(:,2) = VV1(Q+4*I,3) + VV2(Q+4*I,4) + VV5(Q+4*I,4) + VV6(Q+4*I,4) + VV9(Q+4*I,4) + VV10(Q+4*I,4);
end

%VOLTAGE = 0.5
A = mean(VD5(:,1)); % Airspeed Differential Pressure [Pa]
%B = mean(VD5(:,2)); % Aux Differential Pressure [Pa]
%dP = A - B;   %delta P
S_VD5 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(1) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 2
A = mean(V2(:,1));
%B = mean(V2(:,2));
%dP = A - B;   %delta P
S_V2 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(2) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 2.5
A = mean(V2D5(:,1));
%B = mean(V2D5(:,2));
%dP = A - B;   %delta P
S_V2D5 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(3) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 4
A = mean(V4(:,1));
%B = mean(V4(:,2));
%dP = A - B;   %delta P
S_V4 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(4) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 4.5
A = mean(V4D5(:,1));
%B = mean(V4D5(:,2));
%dP = A - B;   %delta P
S_V4D5 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(5) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 6 
A = mean(V6(:,1));
%B = mean(V6(:,2));
%dP = A - B;   %delta P
S_V6 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(6) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 6.5
A = mean(V6D5(:,1));
%B = mean(V6D5(:,2));
%dP = A - B;   %delta P
S_V6D5 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(7) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 8
A = mean(V8(:,1));
%B = mean(V8(:,2));
%dP = A - B;   %delta P
S_V8 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(8) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 8.5
A = mean(V8D5(:,1));
%B = mean(V8D5(:,2));
%dP = A - B;   %delta P
S_V8D5 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(9) = error(A,Tatm,Patm,ratio);

%VOLTAGE = 10
A = mean(V10(:,1));
%B = mean(V10(:,2));
%dP = A - B;   %delta P
S_V10 = VenturiTubeConfig(A, R, Tatm,Patm,ratio);
EVenturiAir(10) = error(A,Tatm,Patm,ratio);

%%
% Airspeed vs. Voltage
VoltageVenturi = [-1 1.5 2 3.5 4 5.5 6 8 9 10];
AirspeedVenturi = [(S_V2-3) S_V2 S_V2D5 S_V4 S_V4D5 S_V6 S_V6D5 S_V8 S_V8D5 S_V10];

figure(1)
errorbar(VoltageVenturi, AirspeedVenturi ,EVenturiAir, 'g');
xlabel('Voltage [V]');
ylabel('Airspeed [m/s]');
title('Average Accuracy of U-tube Manometer vs Pressure Transducer for Both the Pitot Static Probe and Venturi Tube');

save('VenturiAirspeedData.mat', 'VoltageVenturi', 'AirspeedVenturi', 'EVenturiAir');
% 
% plot(VoltageVenturi, AirspeedVenturi);
% xlabel('Voltage [V]');
% ylabel('Airspeed [m/s]');
% title('Venturi Airspeed vs Voltage');

%%
% Obtaining the various pressures with respect to their voltage
for Q = 1:500
	%VOLTAGE = 1.5
    P1D5(:,1) = PP3(Q,3) + PP4(Q,3) + PP7(Q,3) + PP8(Q,3) + PP11(Q,3) + PP12(Q,3);
    %P1D5(:,2) = PP3(Q,4) + PP4(Q,4) + PP7(Q,4) + PP8(Q,4) + PP11(Q,4) + PP12(Q,4);
    %VOLTAGE = 1
    P1(:,1) = PP1(Q,3) + PP2(Q,3) + PP5(Q,3) + PP6(Q,3) + PP9(Q,3) + PP9(Q,3);
    %P1(:,2) = PP1(Q,3) + PP2(Q,4) + PP5(Q,4) + PP6(Q,4) + PP9(Q,4) + PP9(Q,4);
    
    %VOLTAGE = 3.5
    P3D5(:,1) = PP3(Q+I,3) + PP4(Q+I,3) + PP7(Q+I,3) + PP8(Q+I,3) + PP11(Q+I,3) + PP12(Q+I,3);
   % P3D5(:,2) = PP3(Q+I,4) + PP4(Q+I,4) + PP7(Q+I,4) + PP8(Q+I,4) + PP11(Q+I,4) + PP12(Q+I,4);
    %VOLTAGE = 3
    P3(:,1) = PP1(Q+I,3) + PP2(Q+I,3) + PP5(Q+I,3) + PP6(Q+I,3) + PP9(Q+I,3) + PP9(Q+I,3);
    %P3(:,2) = PP1(Q+I,3) + PP2(Q+I,4) + PP5(Q+I,4) + PP6(Q+I,4) + PP9(Q+I,4) + PP9(Q+I,4);
   
    %VOLTAGE = 5.5
    P5D5(:,1) = PP3(Q+2*I,3) + PP4(Q+2*I,3) + PP7(Q+2*I,3) + PP8(Q+2*I,3) + PP11(Q+2*I,3) + PP12(Q+2*I,3);
    %P5D5(:,2) = PP3(Q+2*I,4) + PP4(Q+2*I,4) + PP7(Q+2*I,4) + PP8(Q+2*I,4) + PP11(Q+2*I,4) + PP12(Q+2*I,4);
    %VOLTAGE = 5
    P5(:,1) = PP1(Q+2*I,3) + PP2(Q+2*I,3) + PP5(Q+2*I,3) + PP6(Q+2*I,3) + PP9(Q+2*I,3) + PP9(Q+2*I,3);
    %P5(:,2) = PP1(Q+2*I,3) + PP2(Q+2*I,4) + PP5(Q+2*I,4) + PP6(Q+2*I,4) + PP9(Q+2*I,4) + PP9(Q+2*I,4);

    %VOLTAGE = 7.5
    P7D5(:,1) = PP3(Q+3*I,3) + PP4(Q+3*I,3) + PP7(Q+3*I,3) + PP8(Q+3*I,3) + PP11(Q+3*I,3) + PP12(Q+3*I,3);
    %P7D5(:,2) = PP3(Q+3*I,4) + PP4(Q+3*I,4) + PP7(Q+3*I,4) + PP8(Q+3*I,4) + PP11(Q+3*I,4) + PP12(Q+3*I,4);
    %VOLTAGE = 7
    P7(:,1) = PP1(Q+3*I,3) + PP2(Q+3*I,3) + PP5(Q+3*I,3) + PP6(Q+3*I,3) + PP9(Q+3*I,3) + PP9(Q+3*I,3);
    %P7(:,2) = PP1(Q+3*I,3) + PP2(Q+3*I,4) + PP5(Q+3*I,4) + PP6(Q+3*I,4) + PP9(Q+3*I,4) + PP9(Q+3*I,4);

    %VOLTAGE = 9.5
    P9D5(:,1) = PP3(Q+4*I,3) + PP4(Q+4*I,3) + PP7(Q+4*I,3) + PP8(Q+4*I,3) + PP11(Q+4*I,3) + PP12(Q+4*I,3);
    %P9D5(:,2) = PP3(Q+4*I,4) + PP4(Q+4*I,4) + PP7(Q+4*I,4) + PP8(Q+4*I,4) + PP11(Q+4*I,4) + PP12(Q+4*I,4);
    %VOLTAGE = 9
    P9(:,1) = PP1(Q+4*I,3) + PP2(Q+4*I,3) + PP5(Q+4*I,3) + PP6(Q+4*I,3) + PP9(Q+4*I,3) + PP9(Q+4*I,3);
    %P9(:,2) = PP1(Q+4*I,3) + PP2(Q+4*I,4) + PP5(Q+4*I,4) + PP6(Q+4*I,4) + PP9(Q+4*I,4) + PP9(Q+4*I,4);
end

%VOLTAGE = 1.5
A = mean(P1D5(:,1)); % Airspeed Differential Pressure [Pa]
%B = mean(P1D5(:,2)); % Aux Differential Pressure [Pa]
%dP = A - B;   %delta P
S_P1D5 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(1) = error1(A,Tatm,Patm);

%VOLTAGE = 1
A = mean(P1(:,1));
%B = mean(P1(:,2));
%dP = A - B;   %delta P
S_P1 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(2) = error1(A,Tatm,Patm);

%VOLTAGE = 3.5
A = mean(P3D5(:,1));
%B = mean(P3D5(:,2));
%dP = A - B;   %delta P
S_P3D5 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(3) = error1(A,Tatm,Patm);

%VOLTAGE = 3
A = mean(P3(:,1));
%B = mean(P3(:,2));
%dP = A - B;   %delta P
S_P3 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(4) = error1(A,Tatm,Patm);

%VOLTAGE = 5.5
A = mean(P5D5(:,1));
%B = mean(P5D5(:,2));
%dP = A - B;   %delta P
S_P5D5 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(5) = error1(A,Tatm,Patm);

%VOLTAGE = 5
A = mean(P5(:,1));
%B = mean(P5(:,2));
%dP = A - B;   %delta P
S_P5 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(6) = error1(A,Tatm,Patm);

%VOLTAGE = 7.5
A = mean(P7D5(:,1));
%B = mean(P7D5(:,2));
%dP = A - B;   %delta P
S_P7D5 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(7) = error1(A,Tatm,Patm);

%VOLTAGE = 7
A = mean(P7(:,1));
%B = mean(P7(:,2));
%dP = A - B;   %delta P
S_P7 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(8) = error1(A,Tatm,Patm);

%VOLTAGE = 9.5
A = mean(P9D5(:,1));
%B = mean(P9D5(:,2));
%dP = A - B;   %delta P
S_P9D5 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(9) = error1(A,Tatm,Patm);

%VOLTAGE = 9
A = mean(P9(:,1));
%B = mean(P9(:,2));
%dP = A - B;   %delta P
S_P9 = BernoulliEq(A, R, Tatm,Patm);
EPitotAir(10) = error(A, Tatm,Patm, 1);

%%
% Airspeed vs. Voltage
VoltagePitot = [0 1 2 3 4 5 6 7 8 9];
AirspeedPitot = [S_P1 S_P1D5 S_P3 S_P3D5 S_P5 S_P5D5 S_P7 S_P7D5 S_P9 S_P9D5];

hold on
errorbar(VoltagePitot, AirspeedPitot, EPitotAir, 'm');
xlabel('Voltage [V]');
ylabel('Airspeed [m/s]');
title({'Average Accuracy of U-tube Manometer vs Pressure Transducer', 'for Both the Pitot Static Probe and Venturi Tube'});
legend('U-tube Manometer', 'Pressure Transducer')
set(gca,'Fontsize',20);
save('PitotAirspeedData.mat', 'VoltagePitot', 'AirspeedPitot', 'EPitotAir');

% plot(VoltagePitot, AirspeedPitot);
% xlabel('Voltage [V]');
% ylabel('Airspeed [m/s]');
% title('Pitot Airspeed vs Voltage');

%% Functions
function dv = error(p,Tatm,Patm,ratio)
    R = 8.314;%[J/mol.K]
    delT = 1/4;%[deg C]
    delp = 0.01*6894.76; %[pa]
    delPatm = 3.75;%[1.5% of 250]
    
    dvdp = ((R*Tatm)/(Patm*(1-((ratio)^2))))*sqrt((Patm*(1-((ratio)^2)))/(2*p*R*Tatm));
    dvdTatm = ((p*R)/(Patm*(1-((ratio)^2))))*sqrt((Patm*(1-((ratio)^2)))/(2*p*R*Tatm));
    dvdPatm = (sqrt((Patm*(1-((ratio)^2)))/(2*p*R*Tatm)))*((p*R*Tatm)/(1-((ratio)^2)))/(-(Patm^2));
    dv = sqrt(((dvdp*delp)^2)+((dvdTatm*delT)^2)+((dvdPatm*delPatm)^2));
end

function output = VenturiTubeConfig(deltaP, R, Tatm, Patm, ratio) 
    output = sqrt((2 * deltaP * R * Tatm) / (Patm * (1 - (ratio)^2)));
end


%% Functions
function dv1 = error1(p,Tatm2,Patm2)
    R = 8.314;%[J/mol.K]
    delT = 1/4;%[deg C]
    delp = 0.01*6894.76; %[pa]
    delPatm = 3.75;%[1.5% of 250]
    
    dvdp1 = ((R*Tatm2)/Patm2)/sqrt(2*p*((R*Tatm2)/Patm2));
    dvdTatm1 = ((p*R)/Patm2)/sqrt(2*p*((R*Tatm2)/Patm2));
    dvdPatm1 = (-(p*R*Tatm2)/(Patm2^2))/sqrt(2*p*((R*Tatm2)/Patm2));
    dv1 = sqrt(((dvdp1*delp)^2)+((dvdTatm1*delT)^2)+((dvdPatm1*delPatm)^2));
end
%Function
function output = BernoulliEq(deltaP, R, Tatm, Patm)
     output = sqrt(2*deltaP*((R*Tatm)/Patm));
end 


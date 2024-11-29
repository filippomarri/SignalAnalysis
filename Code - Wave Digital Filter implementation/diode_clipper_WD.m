% Fusari Anna
% Marri Filippo
clear; close all; clc

%% Simscape File for Ground Truth

load('ssc_output.mat')

%% Sampling Frequency
fs = 96e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 1;  % [seconds]

%% Input Signal
% Fundamental Frequency
f0 = 440;

% Time Axis
t = 0:Ts:stop_time;

% Signal Amplitude
A = 1.5;
vin = A * sin(2*pi*f0*t);

% Music Signal Test (Uncomment the following lines for test)
% [vin, fs_rec] = audioread('guitar_input.wav');
% G = 5;                                  % Gain factor 
% vin = G * (vin(:, 1) + vin(:, 2))/2;    % Convert the signal to MONO
% vin = resample(vin, fs, fs_rec);        % Resampling the signal from 44.1 kHz to 96 kHz
% t = 0:Ts:Ts*(length(vin)-1);

%% Circuit Parameters
% Resistive Elements
Rin = 1;
R1 = 1e4;
Rout = 1e4;

% Dynamic Elements
C1 = 1e-6;
C2 = 1e-9;

%% Setting of Free Parameters (Adaptation Conditions)
Z2 = Ts/(2*C2);
Z7 = Rout;
Z6 = Ts/(2*C1);
Z10 = R1;
Z11 = Rin;

%Freed from reflections
Z12 = (Z2 * Z7)/(Z2 + Z7);
Z4 = Z12;
Z9 = Z10 + Z11;
Z8 = Z9;
Z5 = Z6 + Z8;
Z1 = Z5;
Z3 = (Z1 * Z4)/(Z1 + Z4);

%% Computing Scattering Matrices

Bser = [1, 1, 1];
Qpar = [1, 1, 1];
Zser1 = diag([Z9, Z10, Z11]);
Zser2 = diag([Z5, Z6, Z8]);
Zpar1 = diag([Z1, Z3, Z4]);
Zpar2 = diag([Z2, Z7, Z12]);

Sser1= eye(3) - 2*Zser1*Bser'*inv(Bser*Zser1*Bser')*Bser;
Sser2= eye(3) - 2*Zser2*Bser'*inv(Bser*Zser2*Bser')*Bser;
Spar1= 2*Qpar'*inv(Qpar*inv(Zpar1)*Qpar')*Qpar*inv(Zpar1) - eye(3);
Spar2= 2*Qpar'*inv(Qpar*inv(Zpar2)*Qpar')*Qpar*inv(Zpar2) - eye(3);

%% Initialization of Waves
%everything is initialised to zero since both the inductance and the
%capacitance run out of charge and Vin at k = 0 is = 0.
%We highlight that the port are refered to the juunctions
a10 = 0;
b6 = 0;
b2 = 0;
a7 = 0;

%% Initialization of Output Signals
vout = zeros(size(vin));

%% Simulation Algorithm

for n = 1 : length(t)
    %input signal
    a11 = vin(n);
    %elements with memory
    a6 = b6;
    a2 = b2;
    % Forward Scan
    % first layer
    b9 = Sser1(1,:)*[0; a10; a11];
    a8 = b9;
    %second layer
    b5 = Sser2(1,:)*[0; a6; a8];
    a1 = b5;

    b12 = Spar2(3,:)*[a2; a7; 0];
    a4 = b12;
    %third layer
    b3 = Spar1(2,:)*[a1; 0; a4];
    % Local Root Scattering
    % Hint: Use the function 'antiparallel_diodes' to compute the Local Root
    % Scattering
    a3 = antiparallel_diodes(b3, Z3);
    
    % Backward Scan
    %third layer
    b1 = Spar1(1,:)*[a1; a3; a4];
    b4 = Spar1(3,:)*[a1; a3; a4];
    a12 = b4;
    a5 = b1;
    %second layer
    b2 = Spar2(1,:)*[a2; a7; a12];
    b7 = Spar2(2,:)*[a2; a7; a12];

    b6 = Sser2(2,:)*[a5; a6; a8];
    b8 = Sser2(2,:)*[a5; a6; a8];
    a9 = b8;
    %third layer
    b10 = Sser1(2,:)*[a9; a10; a11];
    b11 = Sser1(3,:)*[a9; a10; a11];
    % Read Output
    vout(n) = (a7+b7)/2;

end

% Uncomment the following line to hear the Diode Clipper
 sound(vout, fs)

%% Output Plots

plot_lim = 5/f0; % Limit the plot to just 5 periods of the output signal

figure
set(gcf, 'Color', 'w');
plot(gt(1, :), gt(2, :), 'r', 'Linewidth', 2);
hold on;
plot(t, vout, 'b--', 'Linewidth', 2);
grid on;
xlim([0, plot_lim]);
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [V]','Fontsize',16,'interpreter','latex');
legend('Simscape','WDF','Fontsize',16,'interpreter','latex');
title('Output Signal','Fontsize',18,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, vout - gt(2, :), 'k', 'Linewidth', 2);
grid on;
xlim([0, plot_lim]);
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{out}}$ [V]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',18,'interpreter','latex');

%% Compute Mean Squared Error (MSE)

mse = mean((vout - gt(2, :)).^2);
disp('MSE = ')
disp(mse)

%% Signal fft
VOUT = fft(vout);
VIN = fft(vin);
frequency_axis = 0:fs/length(t): fs - fs/length(t);

figure
set(gcf, 'Color', 'w');
plot(frequency_axis, abs(VIN), 'r', 'Linewidth', 2);
hold on;
plot(frequency_axis, abs(VOUT), 'b', 'Linewidth', 2);
grid on;
xlim([0, 3000]);
xlabel('Frequency [Hz]','Fontsize',16,'interpreter','latex');
ylabel('Abs()','Fontsize',16,'interpreter','latex');
legend('Vin','Vout','Fontsize',16,'interpreter','latex');
title('Fft of the Signals','Fontsize',18,'interpreter','latex');
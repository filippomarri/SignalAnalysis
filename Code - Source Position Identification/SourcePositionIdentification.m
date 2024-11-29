clc;
close all;
clearvars;

%% Load the audio file
[audioData, Fs] = audioread('array_recordings.wav');

L = 0.40; % Total length of the microphones array (in m)

numberOfMicrophones = 16; 

distanceBetweenMic = L / numberOfMicrophones; %Distance between mic

N = length(audioData); %in samples

time_axis = 0 : 1/Fs : (N-1)/Fs;

%Plot the signal captured from each microphone
figure;
for i = 1 : 16
    subplot(4,4,i);
    plot(time_axis, audioData(:,i));
    xlabel('Time (s)');
    title("Microphone "+i);
end

%% Computation of the STFT of audioData

windowLength = 1024;
hopsize = 512; 

%built-in matlab stft
[S_original, f_o, t_o] = stft(audioData, Fs, 'Window', hamming(windowLength), 'OverlapLength', 512);

%our stft
[S, f, t] = ourStft(audioData, Fs, windowLength, hopsize);

%plot both the original and our STFT for comparison
figure;
subplot(2,1,1);
threeDPlot(S(:,:,1), f, t);
title("Our STFT signal")

subplot(2,1,2);
threeDPlot(S_original(:,:,1), f_o, t_o);
title("Built-in matlab STFT signal")

%% Computation of the pseudo-spectrum

soundSpeed = 343; %m/s

%array of angles representing the possible directions of arrival
azimuthAngle = -90:1:90;

pseudoSpectrum = powerEstimation(S, f, t, Fs, numberOfMicrophones, azimuthAngle, distanceBetweenMic, soundSpeed, "geometric", true);

%Computation of the index of the angle related to the maximum value of the
%pseudo-spectrum for each time instant
[~, maxIndex] = max(pseudoSpectrum, [], 2);
maxAngle = azimuthAngle(maxIndex); 

figure;

%pseudo-spectrum plot
subplot(2,1,1);
imagesc(t, azimuthAngle, flipud(pseudoSpectrum.'));
ylabel('DOA [deg]');
xlabel('Time [s]');
title('Pseudo-pectrum');
colorbar;
grid on;

%plot of the values of the angles corresponding to the max power over time
subplot(2,1,2);
plot(t, maxAngle);
xlabel('Time [s]');
ylabel('Max pseudo-spectrum value angle [deg]');
title('Pseudo-spectrum');
grid on;

%% Video

%creation of the video file DOA.avi
v = VideoWriter("DOA.avi");
open(v);

figure;
%plot of the dots representing the microphones array
for i = 0 : numberOfMicrophones - 1
    plot(-L/2 + i*distanceBetweenMic,0,'.b','MarkerSize',25);
    hold on;
end
grid;
xlabel('x [m]');
ylabel('y [m]');

for time = 1 : length(t)

    x_tip = -sind(squeeze(maxAngle(time)));
    y_tip = cosd(squeeze(maxAngle(time)));

    quiver(0,0,x_tip, y_tip);

    frame = getframe(gcf);
    writeVideo(v,frame);
    title("Estimated DOA time: "+ (t(time))+ "[s]");
end

close(v);

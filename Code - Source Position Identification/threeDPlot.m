function threeDPlot(S,f,t)
%The function is used to plot the 3d version of the stft of the signal
%Args:
    %S: the stft of the signal for one of the sixteen channels
    %f: vector representing the frequencies 
    %t: vector representing time instants
    
% Waterfall function plot the S matrix as a third dimension above a grid in
% the (x,y) plane defined by f and t. The edge colours vary according to the
% values specified by S.

    waterfall(f,t,abs(S)'.^2)
    set(gca,XDir="reverse",View=[30 50])
    xlabel("Frequency [Hz]")
    ylabel("Seconds [s]")
end
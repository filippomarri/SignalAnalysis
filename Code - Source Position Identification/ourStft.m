function [S, frequency_axis, time_axis] = ourStft(s, sampling_frequency, window_size, hop_size)
    %The function compute the STFT of the signal 
    %Args:
        % s: input signal
        % sampling_frequency: sampling frequency value
        % window_size: length of the window in samples
        % hop_size: lenght of the non-overlapping section between windows
    %Return:
        %S: the STFT of the signal
        %frequency_axis: array containig the frequency values of the STFT
        %time_axis: array containing the time instants of the STFT
  
    %initialization of an hamming window 
    w = hamming(window_size);

    N = length(s);

    %maximal frame index
    M = floor((N-window_size)/hop_size);
    
    %initialization of the STFT signal
    S = zeros(window_size, M+1, 16);

    %cycling over channels
    for j = 1:16
        s_plane = s(:,j);

        %cycling over the numbers of windows
        for i = 0 : 1 : M

            %cropping the signal and applying the window
            s_cropped = s_plane(i*hop_size + 1 : i*hop_size+window_size).*w;

            %Fourier transform 
            S_cropped = fft(s_cropped);

            %Evaluation of Nyquist frequency in order to centre the STFT
            %in zero
            k = floor((window_size+1)/2);
            X = S_cropped(1:k,:);
            Y = flipud(X);
            S(1:k,i+1,j) = Y;
            S((k+1):2*k,i+1, j) = X;
    
        end
    end

    time_axis = 0 : ((N-1)/sampling_frequency)/M : (N-1)/sampling_frequency;
    time_axis = time_axis.';

    frequency_axis = -(sampling_frequency/2) : sampling_frequency/(2*k) : sampling_frequency/2 - sampling_frequency/(2*k); %up to Nyquist frequency
    frequency_axis = frequency_axis.';
end
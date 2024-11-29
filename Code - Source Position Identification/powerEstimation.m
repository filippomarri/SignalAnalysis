function pseudoSpectrum = powerEstimation(S, f, t, Fs, numberOfMicrophones, azimuthAngle, distanceBetweenMic, soundSpeed, meanComp, normComp)
    %%The function returns the pseudo-spectrum of a sound field given as
    %%input a signal from an array of microphones
 
    %Args:
        %S: STFT of audioData for all 16 channels
        %f: frequncy bins of the STFT
        %t: time instants of the STFT
        %Fs: sampling frequency value
        %numberOfMicrophones: number of microphones (16)
        %azimuthAngle: vector representing the angles from -90° to 90°
        %distanceBetweenMic: distance between microphones in m
        %soundSpeed: the speed of sound
        %meanComp: type of mean computation
            %- geometric
            %- harmonic
            %- arithmetic
        %normComp: boolean value, set it true to normalize values for every time instant

    %Return: 
        %pseudospectrum: a matrix representing the energy content related to
        %a specific direction of arrival (theta from -90° to 90°) and a time
        %instant.
    

    %We consider only the positive frequencies
    k = floor((length(f)+1)/2);
    S_positive = S(1:k,:,:);
    f = f(1:k,:);
    
    %initialization of the pseudo-spectrum
    pseudoSpectrum = zeros(length(t), length(azimuthAngle));

    %cycling over the frequencies
    for freqIndex = 1 : k

        %consider one frequency at a time 
        frequencyAnalysed = f(freqIndex);
        
        %initialization of the power matrix 
        power = zeros(length(t), length(azimuthAngle));
        
        %cycling over the angles
        for a = 1:length(azimuthAngle)

            %consider one DOA at a time
            theta = azimuthAngle(a);
        
            %evaluation of the delay specific for every microphone by
            %considering a certain angle theta
            delay = (0:numberOfMicrophones - 1)*(distanceBetweenMic*sind(theta))*Fs/soundSpeed;
            
            %cycling over time instants
            for i = 1:length(t)
                beamformedSignal = zeros(1, numberOfMicrophones);
                
                %cycling over the number of microphones
                for mic = 1:numberOfMicrophones

                    %saving the STFT of the signal for one channel at a time, considering the current analysed frequency and time instant. 
                    micSignal = S_positive(freqIndex, i, mic);
                  
                    %Appling the delay to the micSignal
                    beamformedSignal(mic) = micSignal * exp(-1i * 2 * pi * frequencyAnalysed * delay(mic) / Fs);
                end
                
                %computing the power of the beamformedSignal for a specific
                %theta and for a certain time instant
                power(i, a) = abs(sum(beamformedSignal))^2;
    
            end
        end
        powerPreAveraged(freqIndex,:,:) = power;
        
    end

    %computation of the mean value of the pseudo-spectrum
    switch meanComp
        case "geometric"
            pseudoSpectrum(:,:) = geomean(powerPreAveraged,1);
        case "harmonic"
            pseudoSpectrum(:,:) = harmmean(powerPreAveraged,1);
        otherwise
            pseudoSpectrum(:,:) = mean(powerPreAveraged,1);
    end

    %computation of the normalization of the pseudo-spectrum
    if normComp == true
        for i = 1 : length(t)
            normVar = max(pseudoSpectrum(i,:));
            pseudoSpectrum(i, :) = pseudoSpectrum(i, :)/normVar;
        end
    end
end

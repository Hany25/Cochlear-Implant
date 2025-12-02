function [x_resampled, fs] = MTE252P1(audioFile)
% Phase 1
% Input: audio file filename
% Output: plots (input sound and output cos) and sound files

% Read Audio File
[x, fs] = audioread(audioFile); 

% Check if converting to mono needed
if size(x,2) == 2
    x = mean(x, 2); % Average the two channels
else
    disp('Mono signal detected.');
end

% Playing audio
disp('Playing the input sound...');
% sound(x, fs);
% pause(length(x)/fs + 1); % Wait until playback ends

% Write a new file for the output sound
newFile = 'output_mono.wav';
audiowrite(newFile, x, fs);

% Downsampling 
target_fs = 16000;
if fs ~= target_fs
    x_resampled = resampleAudio(x, fs, target_fs);
    fs = target_fs;
else
    x_resampled = x;
end

function y = resampleAudio(x, fs_in, fs_out)
    t_in = (0:length(x)-1) / fs_in;
    t_out = 0:1/fs_out:(length(x)-1)/fs_in;
    y = interp1(t_in, x, t_out, 'linear')';
end

% Generate a cosine signal 
N = length(x_resampled);
t = (0:N-1)/fs; % Time vector
f0 = 1000; % 1kHz tone
cos_signal = cos(2*pi*f0*t)';

% Play the signal
disp('Playing 1 kHz cosine signal...');
%sound(cos_signal, fs);
%pause(length(cos_signal)/fs + 1);

% Plot the audio 
figure;
plot(x);
xlabel('Sample Number'); ylabel('Amplitude');
title('Waveform of Input Sound');

% % Plot two cycles of cosine
% T = 1/f0; % Period
% samples_per_cycle = round(T * fs);
% figure;
% plot(t(1:2*samples_per_cycle), cos_signal(1:2*samples_per_cycle));
% xlabel('Time (s)'); ylabel('Amplitude');
% title('Cos Wave 1 kHz for Two Cycles');
% 
% % Save the cosine tone
% audiowrite('cos_1kHz.wav', cos_signal, fs);

end
function FilterBankPhase2(x, fs)
% Phase 2
% Input: x - preprocessed, resampled mono audio
%        fs - sampling rate (e.g., 16 kHz)

% Use via:
% [x, fs] = MTE252P1('recording.wav');
% MTE252P2(x, fs);

% Filter bank design
N = 8;  % Number of channels
f_low = 100;
f_high = 8000;
band_edges = logspace(log10(f_low), log10(f_high), N+1); % logarithmic spacing

filters = cell(N,1);
for i = 1:N
    f1 = band_edges(i);
    f2 = band_edges(i+1);
    
    % Design bandpass filter (4th-order Butterworth)
    [b, a] = butter(4, [f1 f2]/(fs/2), 'bandpass');
    filters{i} = struct('b', b, 'a', a, 'range', [f1 f2]); % Stores coefficients and bandrange
end

% Filter the signal
filtered_signals = cell(N,1);
for i = 1:N
    filtered_signals{i} = filter(filters{i}.b, filters{i}.a, x);
end

% Plot lowest and highest frequency channels
figure;
subplot(2,1,1);
plot(filtered_signals{1});
title(sprintf('Lowest Band (%.0f–%.0f Hz)', filters{1}.range));
xlabel('Samples'); ylabel('Amplitude');

subplot(2,1,2);
plot(filtered_signals{N});
title(sprintf('Highest Band (%.0f–%.0f Hz)', filters{N}.range));
xlabel('Samples'); ylabel('Amplitude');

% Rectification
rectified_signals = cellfun(@abs, filtered_signals, 'UniformOutput', false); % abs value rectification sample-wise

% Envelope extraction
fc = 400; % cutoff freq for lowpass
[b_lp, a_lp] = butter(4, fc/(fs/2), 'low');

envelopes = cell(N,1);
for i = 1:N
    envelopes{i} = filter(b_lp, a_lp, rectified_signals{i});
end

% Plot envelopes of lowest and highest bands
figure;
subplot(2,1,1);
plot(envelopes{1});
title(sprintf('Envelope - Lowest Band (%.0f–%.0f Hz)', filters{1}.range));
xlabel('Samples'); ylabel('Amplitude');

subplot(2,1,2);
plot(envelopes{N});
title(sprintf('Envelope - Highest Band (%.0f–%.0f Hz)', filters{N}.range));
xlabel('Samples'); ylabel('Amplitude');

disp('Phase 2 completed successfully.');

end
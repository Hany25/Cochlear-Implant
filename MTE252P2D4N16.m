function MTE252P2D4N16(x, fs)
% Phase 2
% Input: x - preprocessed, resampled mono audio
%        fs - sampling rate (e.g., 16 kHz)

% Use via:
% [x, fs] = MTE252P1('recording-24s.wav');
% MTE252P2D4N16(x, fs);

% Filter bank design parameters
N = 16;  % Number of channels
f_low = 100; % Hz
f_high = 7999; % less than 8000 to stay within Nyquist bounds
band_edges = logspace(log10(f_low), log10(f_high), N+1); % logarithmic spacing
filter_order = 4;

filters = cell(N,1);
for i = 1:N
    f1 = band_edges(i);
    f2 = band_edges(i+1);

    % Bandpass filter - 4th-order Butterworth
    [b, a] = butter(filter_order, [f1 f2]/(fs/2) , 'bandpass');
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
fc = 800; % cutoff freq for lowpass
norm_cutoff = fc / (fs/2); % normalize to Nyquist frequency

[b_lp, a_lp] = butter(filter_order, norm_cutoff, 'low');

envelopes = cell(N,1);
for i = 1:N
    % recursive difference equation implementation
    x_in = rectified_signals{i};
    y_out = zeros(size(x_in));
    for n = filter_order+1:length(x_in)
        y_out(n) = b_lp(1)*x_in(n) ...
                 + b_lp(2)*x_in(n-1) ...
                 + b_lp(3)*x_in(n-2) ...
                 + b_lp(4)*x_in(n-3) ...
                 + b_lp(5)*x_in(n-4) ...
                 - a_lp(2)*y_out(n-1) ...
                 - a_lp(3)*y_out(n-2) ...
                 - a_lp(4)*y_out(n-3) ...
                 - a_lp(5)*y_out(n-4);
    end
    envelopes{i} = y_out;
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

% PHASE 3 - Fc = 800 Hz
% Determine center frequencies
centers = zeros(N,1);
for i = 1:N
    centers(i) = sqrt(filters{i}.range(1) * filters{i}.range(2)); % Geometric mean since log spacing
end

L = length(x);
t = (0:L-1)'/fs; % make column vector for ease of use later

synth_channels = cell(N,1);

for i = 1:N
    carrier = cos(2*pi*centers(i)*t);   % cosine at channel freq
    env = envelopes{i};
    synth_channels{i} = carrier .* env; % amplitude modulation
end

% Add all the channels together
y = zeros(L,1);
for i = 1:N
    y = y + synth_channels{i};
end

% Normalize signal
maxAmp = max(abs(y));
y = y / maxAmp;

% Bonus task: quantitative comparison
band_error = zeros(N,1);
filtered_y = cell(N,1);
for i = 1:N
    filtered_y{i} = filter(filters{i}.b, filters{i}.a, y);

    Ein  = sqrt(mean(filtered_signals{i}.^2)); % input band RMS
    Eout = sqrt(mean(filtered_y{i}.^2)); % output band RMS

    if Ein > 0
        band_error(i) = abs(Eout - Ein) / Ein;
    else
        band_error(i) = 0;
    end
end

overall_error = mean(band_error);
similarity_pct = max(0, 1 - overall_error) * 100; % Outputs a percentage
fprintf('Bonus Task: Similarity = %.1f%%\n', similarity_pct);

% Play and save output
disp("Playing synthesized Phase 3 output...");
sound(y, fs);
pause(length(y)/fs + 0.5);

audiowrite("Phase3_Output.wav", y, fs);
disp("Saved as Phase3_Output.wav");

% Plot comparison
figure;
subplot(2,1,1); plot(x); title("Input Waveform");
subplot(2,1,2); plot(y); title("Phase 3 Output Waveform");

end
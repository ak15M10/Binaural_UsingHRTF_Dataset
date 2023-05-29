
% The zero-crossing rate is simply the
% number of times the signal crosses the
% zero-axis per second, 
% divided by the length of the signal.
Fs=fs01;
t=0:1/Fs:5;
music11=2*cos(2*pi*500*t+pi/3);

zc_rate = mean(abs(diff(sign(music11))))/(2*fs01);


% Compute envelope of audio signal
audio_env = abs(hilbert(music11));

% Plot audio signal and its envelope
t = (0:length(music11)-1)/Fs;
figure(1)
subplot(2,1,1);
plot(t, music11);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Audio Signal');

subplot(2,1,2);
plot(t, audio_env);
xlabel('Time (s)');
ylabel('Envelope Amplitude');
title('Envelope of Audio Signal');
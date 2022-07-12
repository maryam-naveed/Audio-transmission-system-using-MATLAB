%%%% Reading and Plotting Audio Signal with Noise %%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; close all
[signal_orignal_1,Fs] = audioread('test.wav');
[signal_orignal_2,Fs] = audioread('test1.wav');
samples = 2*Fs; % Fs = 48000
T = 1/Fs;
N = samples; % N = number of samples
n = 0:samples-1; % n = 0:N-1
t = n*T;  % t = nT = n/Fs
noisy_signal_1 = signal_orignal_1(1:samples); % sampling signal 1
noisy_signal_2 = signal_orignal_2(1:samples); % sampling signal 2
f = -Fs/2:Fs/samples:Fs/2-(Fs/samples); % initializing f

% Noisy signal after lowpass filtering

f_cut = 2000; % filter cutoff frequency in Hz
wn = f_cut/(Fs/2);
[num,den] = butter(10, wn); % butterworth filter design
filtered_noisy_signal_1 = filter(num,den,noisy_signal_1); % filters the input data using a rational
%     transfer function defined by the numerator and denominator
%     coefficients num and den
filtered_noisy_signal_2 = filter(num,den,noisy_signal_2);
sound(filtered_noisy_signal_1,Fs) % plays data as sound
pause(3) % stops execution temporarily
sound(filtered_noisy_signal_2,Fs) % plays data as sound
m1 = filtered_noisy_signal_1; % filtered signals are the message signals
m2 = filtered_noisy_signal_2;

% plotting message signals (time domain and frequency domain)

figure % creates figure window
subplot(221) % subplot(row,column,index) - 2 rows, 2 columns, first index
plot(t,m1) % plotting message signal m1
title('Time Domain Message Signal (m1)'); % adds title
grid on % turns on the grid
subplot(222)
plot(f,abs(fftshift(fft(m1)))) % FFT (not normalized)
title('Magnitude Spectrum of Noisy Signal (m1)');
grid on
xlim([-3000 3000]) % sets x-axis limits
subplot(223)
plot(t,m2) % plotting message signal m2
title('Time Domain Message Signal (m2)');
grid on
subplot(224)
plot(f,abs(fftshift(fft(m2))))    % FFT (not normalized)
title('Magnitude Spectrum of Noisy Signal (m2)');
xlim([-3000 3000])
grid on

fc1 = 7500; % 5kHz + 2500
c1 = cos(2*pi*fc1*t); % initializing first carrier signal 

% plotting first carrier signal (frequency domain and time domain)

figure(2) % creates figure 2 window
subplot(221)
plot(t,c1,'linewidth',2); % plots c(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('|c1(t)|','fontweight','bold'); % adds label along y-axis
title('c1(t)') % adds title
xlim([0 0.005]); % sets x-axis limits
ylim([-1.5 1.5]);

Xk = fft(c1)/samples; % computes fft
Y = abs(Xk); % computes the magnitude
X = fftshift(abs(Xk)); % shifts to center of spectrum
fx = -Fs/2:Fs/length(X):Fs/2-(Fs/length(X));
subplot(222)
stem(fx,X,'filled','linewidth',2); % plots discrete sequence data
grid on;
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X2(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency domain spectrum of c1(t)'); % adds title
ylim([0 1]);

fc2 = 12500; % 7500 + 5000
c2 = cos(2*pi*fc2*t); % initializing first carrier signal

% plotting second carrier signal (frequency domain and time domain)

subplot(223)
plot(t,c2,'linewidth',2); % plots c(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('|c2(t)|','fontweight','bold'); % adds label along y-axis
title('c2(t)') % adds title
xlim([0 0.005]);
ylim([-1.5 1.5]);

Xk1 = fft(c2)/samples; % computes fft
Y1 = abs(Xk1); % computes the magnitude
X1 = fftshift(abs(Xk1)); % shifts to center of spectrum
fx1 = -Fs/2:Fs/length(X1):Fs/2-(Fs/length(X1));
subplot(224) % plots at the 4th index (2 rows, 2 columns)
stem(fx1,X1,'filled','linewidth',2);
grid on;
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X2(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency domain spectrum of c1(t)'); % adds title
ylim([0 1]);

% modulated signal
% performing SSB modulation using message signal 1 and carrier 1
u1 = m1'.*cos(2*pi*fc1*t)+imag(hilbert(m1')).*sin(2*pi*fc1*t); 

% plotting first modulated signal (frequency domain and time domain)

figure(3) % creates figure 3 window
subplot(221) % plots the graph at the specified position
plot(t,u1,'linewidth',2); % plots m(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('|u1(t)|','fontweight','bold'); % adds label along y-axis
title('u1(t) = m1(t)*c1(t)') % adds title
Xk2 = fft(u1)/N; % computes fft
Y2 = abs(Xk2); % computes the magnitude
X2 = fftshift(abs(Xk2)); % shifts to center of spectrum
fx2 = -Fs/2:Fs/length(X2):Fs/2-(Fs/length(X2));
subplot(222)
plot(fx2,X2);
grid on;
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X2(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency Domain Spectrum of u1'); % adds title

% performing SSB modulation using message signal 2 and carrier 2

u2 = m2'.*cos(2*pi*fc2*t)+imag(hilbert(m2')).*sin(2*pi*fc2*t); 

% plotting second modulated signal (frequency domain and time domain)

subplot(223) % plots the graph at the specified position
plot(t,u2,'linewidth',2); % plots m(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('|u2(t)|','fontweight','bold'); % adds label along y-axis
title('u2(t) = m2(t)*c2(t)') % adds title
Xk3 = fft(u2)/N; % computes fft
Y3 = abs(Xk3); % computes the magnitude
X3 = fftshift(abs(Xk3)); % shifts to center of spectrum
fx3 = -Fs/2:Fs/length(X3):Fs/2-(Fs/length(X3));
subplot(224)
plot(fx3,X3);
grid on;
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X3(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency Domain Spectrum of u2'); % adds title

% multiplexing both modulated message signals

u = u1 + u2; 

% plotting the multiplexed signal (time and frequency domain)

figure(4) % craetes figure 4 window
subplot(211) % plots the graph at the specified position
plot(t,u,'linewidth',2); % plots u(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('|u(t)|','fontweight','bold'); % adds label along y-axis
title('u(t) = u1+u2') % adds title
Xk4 = fft(u)/N; % computes fft
Y4 = abs(Xk4); % computes the magnitude
X4 = fftshift(abs(Xk4)); % shifts to center of spectrum
fx4 = -Fs/2:Fs/length(X4):Fs/2-(Fs/length(X4));
subplot(212)
plot(fx4,X4);
grid on;
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X4(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency Domain Spectrum of u(t)'); % adds title

% demodulation of received signal for the first message signal

r1 = u.*c1; % demodulation using c1
wn = 2500/(Fs/2); % wn = cut-off frequency/(Fs/2)
[a,b] = butter(10,wn); % designs butterworth filter
mest_1 = filter(a,b,r1); % filers the function

% plotting the demodulated signal 1 (time and frequency domain)

figure(5)
subplot(421) % plots the graph at the specified position
plot(t,r1,'linewidth',2); % plots m(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('|r1(t)|','fontweight','bold'); % adds label along y-axis
title('r1(t) = u(t)*c1(t)') % adds title
Xk5 = fft(r1)/N; % computes fft
Y5 = abs(Xk5); % computes the magnitude
X5 = fftshift(abs(Xk5)); % shifts to centre of spectrum
fx5 = -Fs/2:Fs/length(X5):Fs/2-(Fs/length(X5));
subplot(422)
plot(fx5,X5);
grid on; % turns on the grid
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X5(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency Domain Spectrum of r1(t)'); % adds title
subplot(423) % plots the graph at the specified position
plot(t,mest_1,'linewidth',2); % plots m(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('mest 1','fontweight','bold'); % adds label along y-axis
title('DC Block');
Xk6 = fft(mest_1)/N; % computes fft
Y6 = abs(Xk6); % computes the magnitude
X6 = fftshift(abs(Xk6)); % shifts to centre of spectrum
fx6 = -Fs/2:Fs/length(X6):Fs/2-(Fs/length(X6));
subplot(424) % plots the graph at the sepcified position
plot(fx6,X6);
grid on;
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X6(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency Domain Spectrum (DC)'); % adds title
xlim([-3000 3000])

% demodulation of received signal for the second message signal

r2 = u.*c2; % demodulation using c2
wn = 2500/(Fs/2); % wn = cut-off frequency/(Fs/2)
[a,b] = butter(10,wn); % designs butterworth filter
mest_2 = filter(a,b,r2); % filters the data using rational transfer function
figure(5)
subplot(425) % plots the graph at the specified position
plot(t,r2,'linewidth',2); % plots m(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('|r2(t)|','fontweight','bold'); % adds label along y-axis
title('r2(t) = u(t)*c2(t)') % adds title
Xk7 = fft(r2)/N; % computes fft
Y7 = abs(Xk7); % computes the magnitude
X7 = fftshift(abs(Xk7)); % shifts to centre of spectrum
fx7 = -Fs/2:Fs/length(X7):Fs/2-(Fs/length(X7));
subplot(426)
plot(fx7,X7);
grid on;
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X7(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency Domain Spectrum of r2(t)'); % adds title
subplot(427) % plots the graph at the specified position
plot(t,mest_2,'linewidth',2); % plots m(t)
grid on; % turns on the grid
xlabel('t','fontweight','bold'); % adds label along x-axis
ylabel('mest 2','fontweight','bold'); % adds label along y-axis
title('DC Block');
Xk8 = fft(mest_2)/N; % computes fft
Y8 = abs(Xk8); % computes the magnitude
X8 = fftshift(abs(Xk8)); % shifts to centre of spectrum
fx8 = -Fs/2:Fs/length(X8):Fs/2-(Fs/length(X8));
subplot(428)
plot(fx6,X6);
grid on; % turns on the grid
xlabel('-Fs/2 to Fs/2','fontweight','bold'); % adds label along x-axis
ylabel('|X8(k)|','fontweight','bold'); % adds label along y-axis
title('Frequency Domain Spectrum (DC)'); % adds title
xlim([-3000 3000]) % sets x-axis limits

sound(mest_1,Fs) % plays data as sound
pause(3) % stops execution temporarily
sound(mest_2,Fs) % plays data as sound

close all; clear all
[signal_orignal_1,Fs] = audioread('test1.wav'); % reads the audio file
msg = signal_orignal_1;
samples = 2*Fs; % Fs = 48000
T = 1/Fs;
n = 0:samples-1; % n = 0:N-1
t = n*T;  % t = nT = n/Fs
noisy_signal_1 = signal_orignal_1(1:samples); % sampling signal 1
f_cut = 2500; % filter cutoff frequency in Hz
wn = f_cut/(Fs/2);
[num,den] = butter(10, wn); % butterworth filter design
filtered_noisy_signal_1 = filter(num,den,noisy_signal_1);
sound(filtered_noisy_signal_1,Fs) % plays data as sound

v = 8;
N = 2^v; % quantization levels
x_max = max(abs(filtered_noisy_signal_1)); % max amp
x_min = min(abs(filtered_noisy_signal_1));
del = (2*x_max)/N; % computing del
partition = [x_min:del:x_max]; % partitions data
codebook = [x_min-del/2:del:x_max+del/2];
[ind,quants] = quantiz(filtered_noisy_signal_1,partition,codebook); % returns a vector that tells which interval each input is in
pause(3) % adds a pause of 3 seconds
sound(quants,Fs) % plays data as sound

if(max(ind) > N-1)
    for i = 1:length(ind)
    ind(i) = ind(i) - 1;
    end
end
for i = 1:length(quants)
if(quants(i) == x_max-del/2)
    quants(i) = x_min+del/2;
end
end

ind = abs(ind);

subplot(311); % subplot(row,column,index) - 3 rows, 1 column, first index
plot(t,msg,'linewidth',2); % plots the message signal wrt t
title('Message Signal'); % adds title
xlabel('t','fontweight','bold')
ylabel('amp','fontweight','bold')
grid on; % turns on the grid
subplot(312);
plot(t,quants,'linewidth',2);
title('Quantized Signal');
xlabel('t','fontweight','bold')
ylabel('amp','fontweight','bold')
grid on;
xlim([0 0.1])

ybin = de2bi(ind,'left-msb');
vec = reshape(ybin',1,[]); % converts code matrix to vector
subplot(313);
stairs(vec,'linewidth',2); % creates a stairstep graph
title('PCM PLOT');
xlabel('t','fontweight','bold')
ylabel('amp','fontweight','bold')
xlim([18000 20000]);
grid on

% BPSK
data_bits = vec;
N1 = 3;  % samples per symbol
Tb = 1;  % Bit duration
f = 50000000; % f in Hz
Ab = 1;  %Bit Amplitude
t1 = 0:1/N1:Tb-(1/N1);
w = sqrt(2/Tb)*cos(2*pi*f*t1); % Defining the Basis Function
pulse_bit_1 = Ab*w;
pulse_bit_0 = -1*Ab*w;
BPSK_Pulse = [];
a = length(data_bits)

for i = 1:length(data_bits);
     if data_bits(i) == 1;
         BPSK_Pulse = [BPSK_Pulse pulse_bit_1];
     else
         BPSK_Pulse = [BPSK_Pulse pulse_bit_0];
         a = i
     end
end

figure(2)
plot([0:length(BPSK_Pulse)-1]/N1,BPSK_Pulse,'linewidth',2)
xlabel('time','fontweight','bold')
ylabel('amplitude','fontweight','bold')
title('BPSK','fontweight','bold') % adds title
grid on
xlim([10010 10090])

variance = 1;
r=BPSK_Pulse+sqrt(variance)*randn(1,length(BPSK_Pulse)); % Addition of Noise in the Channel
% Receiver (Bits Detection using Correlator)
Rcvd_Bits = [];
Eb =sum(pulse_bit_1.^2); %Bit Energy
Ew = sum(w.^2);  % Energy of Basis Function

for i=1:N1:length(BPSK_Pulse)
z = 1/sqrt(Eb*Ew)*sum(w.*r(i:i+(N1-1)));  % Computing Correlation
% Multiplying with 1/sqrt(Eb*Ew) to normalize the correlator output.
if z > 0
Rcvd_Bits = [Rcvd_Bits 1];
else
Rcvd_Bits = [Rcvd_Bits 0];
end
end

vec_re = reshape(Rcvd_Bits',v,[])';
index = bi2de(vec_re,'left-msb'); % converts binary numbers to decimal numbers
q_rec = del*index'+x_min+del/2; % getting back the quantized values
figure(3)
plot(t,q_rec,'linewidth',2)
title('Demodulated Signal')
xlabel('t','fontweight','bold')
ylabel('amp','fontweight','bold')
grid on;
xlim([0 0.1])
sound(q_rec,Fs)
clear;
close all

load('data/CSI_breathe.mat');

% CSI_sanit_pha = CSI_sanit_pha(1:7500);

figure('Position',[200 600 1500 300])
plot(CSI_sanit_pha)
ylabel('phase')

% figure(2)
% for i = 21:length(CSI_sanit_pha)
%     i
%     
%     scatter(real(CSI_sanit_vec(i)),imag(CSI_sanit_vec(i)),'filled','b');    
%     xlim([-0.1 0.1]);ylim([-0.1 0.1])
%     hold on
%     
%     waitforbuttonpress();
% end

%%
T_sample = 5e-3;
Fs = 1/T_sample;

% figure(1)
% highpass(CSI_sanit_pha,0.33,Fs);
CSI_sanit_pha = highpass(CSI_sanit_pha,0.33,Fs);

% figure(2)
CSI_sanit_pha = lowpass(CSI_sanit_pha,0.67,Fs);

figure('Position',[200 600 1500 300])
plot(CSI_sanit_pha)
ylabel('phase')

[f,s] = get_fft(CSI_sanit_pha,Fs);

[max_s,I_s] = max(abs(s));

breath_freq = f(I_s)
breath_oneMin = breath_freq*60

figure('Position',[200 200 1500 300])
stem(f,abs(s))
xlim([0 10])

%% FFT
function [f,P1] = get_fft(x, Fs)
    Y = fft(x);
    L = length(x);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))/L;
end


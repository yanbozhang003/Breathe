clear;
close all

load('../csi_log/yanbo/test/breathe.mat');

num_CSI = length(CSI_st_all);

CSI_sanit_vec = zeros(1,1);
for i = 1:1:num_CSI
    
    CSI = CSI_st_all(i).csi;
    
    csi_11 = mean(CSI(1,1,:)./CSI(1,2,:));
    
    CSI_sanit_vec(i,1) = csi_11;
end

CSI_sanit_pha = angle(CSI_sanit_vec);

figure(1)
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

[f,s] = get_fft(CSI_sanit_pha,Fs);

s(1) = 0;
[max_s,I_s] = max(abs(s));

breath_freq = f(I_s)
breath_oneMin = breath_freq*60

figure(2)
stem(f,abs(s))


%% FFT
function [f,P1] = get_fft(x, Fs)
    Y = fft(x);
    L = length(x);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))/L;
end


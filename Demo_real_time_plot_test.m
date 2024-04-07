function collect_csi

clear all;

cleanupObj = onCleanup(@cleanMeUp);

addpath('./msocket_func');

server_sock 	= mslisten(6767);
sock_vec_in     = zeros(1,1,'int32');

sock_vec_in(1,1)= server_sock;

sock_cnt_in    = length(sock_vec_in);
sock_max    = max(sock_vec_in);
sock_min    = min(sock_vec_in);
car_num = 56;

total = 0;
csi_cnt = 0;

PLOT = 0;
n_sample_file = 5000000;

T_sample = 5e-2;
Fs = 1/T_sample;

CSI_buffer = zeros(1,1);
T_window = 10;       % window length is 5 sec.
L_window = T_window*Fs;
t_sample = 0:1/Fs:L_window/Fs-1/Fs;

breathe_est_interval = 20;      % breathe estimation interval is 10 sec.
L_est_interval = breathe_est_interval*Fs;
bpf_low = 0.025;
bpf_high = 0.067;
[b,a] = butter(6,[bpf_low,bpf_high]/(Fs/2)); % 6-order butterworth filter
t_est_sample = 0:1/Fs:L_est_interval/Fs-1/Fs;

figure('Position',[200 600 1500 800])

while 1
    tic
    
    [sock_vec_out,sock_cnt_out,CSI_struct] = msCSI_server_tmp(sock_cnt_in,sock_vec_in,sock_min,sock_max,2);
    if (length(sock_vec_in) > 1)                        % we connect with at least 1 client
        if (length(CSI_struct) >= 1)                     % the output CSI structure must not be empty
            for csi_st_idx = 1:1:length(CSI_struct)     % all the structure
                CSI_entry   = CSI_struct(1,csi_st_idx);
                N_tx        = CSI_entry.nc;
                N_rx        = CSI_entry.nr;
                num_tones   = CSI_entry.num_tones;
                pay_len     = CSI_entry.payload_len;
%               
                if N_rx < 3  || num_tones~= car_num
                    continue;
                end
              
                if N_tx < 3  || num_tones~= car_num
                    continue;
                end
              
                if CSI_struct.MAC_idx ~= 14562
                    continue;
                end
                
                if CSI_struct.noise ~= 0
                    continue;
                end

                if isempty(CSI_struct.csi) == 1
                    continue;
                end
                
                total = total + 1

                %% csi sanit
                csi = CSI_struct.csi;

                csi_11 = mean(csi(1,1,:)./csi(1,2,:));
                csi_33 = mean(csi(3,3,:)./csi(3,2,:));

                %% breathe estimation
                CSI_buffer(total,1) = csi_11;

                if total >= L_window
                    length_buffer = length(CSI_buffer);
                    csi_plot = CSI_buffer(length_buffer-L_window+1:end);
                    pha_plot = unwrap(angle(csi_plot));

                    subplot(3,1,1)
                    plot(t_sample,pha_plot);
%                     ylim([0 pi]);
                    xlabel('Time (s)'); ylabel('Phase (unwrapped)');
                    drawnow()
                    
                    if mod(total,L_est_interval) == 0
                        [freq,spec,breath_oneMin,pha_filter]=breath_est(unwrap(angle(CSI_buffer(total-L_est_interval+1:end))),Fs,b,a);
                        
                        subplot(3,1,2)
                        plot(t_est_sample,pha_filter)
                        xlabel('Time (s)'); ylabel('Phase (unwrapped)');
                        drawnow()
                        
                        subplot(3,1,3)
                        stem(freq,abs(spec))
                        xlim([0 10])
                        xlabel('Frequency (Hz)'); ylabel('Amplitude');
                        title(['breathe: ', num2str(breath_oneMin),' times/min']);
                        drawnow()
                    end
                end
                
                if mod(total, n_sample_file) == 0 
                    cleanMeUp();
                    break;
                end
            end
        end
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%            adjust the sockets 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sock_cnt_in     = sock_cnt_out;
    sock_vec_in     = sock_vec_out;
    sock_max        = max(sock_vec_in);
    sock_min        = min(sock_vec_in);    
    
    toc
end

%% FFT
function [f,P1] = get_fft(x, Fs)
    Y = fft(x);
    L = length(x);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))/L;
end

%% breathe estimation
function [f,s,breath_oneMin,pha_filter] = breath_est(pha_sanit,Fs,b,a)
%     pha_filter= filter(b,a,pha_sanit);
    pha_filter = highpass(pha_sanit,0.25,Fs);
    pha_filter = lowpass(pha_filter,0.67,Fs);

    [f,s] = get_fft(pha_filter,Fs);
    
    [max_s,I_s] = max(abs(s));

    breath_freq = f(I_s);
    breath_oneMin = breath_freq*60;
end

%% clean up socket
function cleanMeUp()
    for i = 1:1:length(sock_vec_out)
        fprintf('close all active socket and exit!!\n');
        msclose(sock_vec_out(i,1));
    end
    close all;
    clear all;
end

end



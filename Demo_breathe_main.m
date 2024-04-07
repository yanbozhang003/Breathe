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

index = 0;
total = 0;
FILL_FLAG = 0;
csi_cnt = 0;
edge_count = 0;

PLOT = 0;
n_sample_file = 300000;

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
                
                if pay_len ~= 120
                    edge_count = edge_count+1;
                    if mod(edge_count,2)==1
                        FILL_FLAG = 1;
                    elseif mod(edge_count,2)==0
                        FILL_FLAG = 0;
                    end
                end
                
                if FILL_FLAG == 1
                    close all
                    % start filling CSI_st_all
                    csi_cnt = csi_cnt + 1;
                    CSI_st_all(csi_cnt) = CSI_entry;
                elseif FILL_FLAG == 0
                    if csi_cnt ~= 0
                     %% csi sanit
                        csi_cnt  = 0;
                        csi_11 = zeros(1,1); csi_33 = zeros(1,1); edge = zeros(1,1);
                        edge_cnt = 2; edge(1) = 1; edge(2) = length(CSI_st_all);
                        for pkt_idx = 1:1:length(CSI_st_all)
                            csi = CSI_st_all(1,pkt_idx).csi;
                            
                            csi_cnt = csi_cnt + 1
                            csi_11(csi_cnt) = mean(csi(1,1,:)./csi(1,2,:));
                            csi_33(csi_cnt) = mean(csi(3,3,:)./csi(3,2,:));
                        end
                     
                     %% breathe estimation
                        CSI_sanit_pha = unwrap(angle(csi_11));
                        
                        save('data/CSI_breathe.mat','CSI_sanit_pha');
                        %% compute frequency
                        T_sample = 5e-3;
                        Fs = 1/T_sample;
                        t_sample = 0:1/Fs:length(CSI_sanit_pha)/Fs-1/Fs;
                        
                        figure('Position',[200 600 1500 900])
                        subplot(3,1,1)
                        plot(t_sample,CSI_sanit_pha)
                        xlabel('Time (s)'); ylabel('Phase (unwrapped)');
%                         drawnow()
                        
                        CSI_sanit_pha = highpass(CSI_sanit_pha,0.25,Fs);
                        CSI_sanit_pha = lowpass(CSI_sanit_pha,0.67,Fs);
                        
%                         figure('Position',[200 600 1500 300])
                        subplot(3,1,2)
                        plot(t_sample,CSI_sanit_pha)
                        xlabel('Time (s)'); ylabel('Phase (unwrapped)');
%                         drawnow()

                        [f,s] = get_fft(CSI_sanit_pha,Fs);

                        [max_s,I_s] = max(abs(s));

                        breath_freq = f(I_s)
                        breath_oneMin = breath_freq*60

%                         figure('Position',[200 200 1500 300])
                        subplot(3,1,3)
                        stem(f,abs(s))
                        xlim([0 10])
                        xlabel('Frequency (Hz)'); ylabel('Amplitude');
                        title(['breathe: ', num2str(breath_oneMin),' times/min']);
                        drawnow()
                        
                        % clear CSI_st_all
                        clear CSI_st_all;
                        csi_cnt = 0;
                    end
                end
                
                total = total + 1
                
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

    
function cleanMeUp()
    for i = 1:1:length(sock_vec_out)
        fprintf('close all active socket and exit!!\n');
        msclose(sock_vec_out(i,1));
    end
    close all;
    clear all;
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

end



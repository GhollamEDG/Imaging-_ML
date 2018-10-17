% multi_AoA_3D
close all;
clc;    clear all;
addpath('../function');

%% Processing selection
IF = struct('debug',0, ...
            'DC_blk',0, ...
            'normal',1, ...
            'save_bg',0, ...
            'bg_sub',1, ...
            'slow_ramp',0, ...
            'cut_lh_freq',1, ...
            'save',0, ...
            'load',0);

%% Parameters
scan_mode = 'cross';
% scan_mode = 'rect';
wvform = 'fsr_tri_6ms';
ant_ref = 'o';

date_ref = '1009';      % reference of experiment date
obj_ref = '23';
run_ref = 1;

Rs = 5e5;   % USRP sampling rate 
T_mov = 30;
T_meas = 35;            % time of sampling in seconds
v_ls = 1.25e-2;
N_spl = T_meas * Rs;    % number of samples

var_lib_3D;

int_cbl_len = 2.14;

dist_ulim = 10;
clim = 100000000;

file_addr = sprintf('~/Jayden/Imaging_ML/data/%s/o%s_r%d',date_ref,obj_ref,run_ref);
% file_addr = sprintf('~/Dropbox/Project/Imaging with ML/data/%s/o%s_r%d',date_ref,obj_ref,run_ref);

% Cut off USRP setup
RX_start = Rs/10;    % skip the first half second samples due to USRP start up
N_spl = N_spl-RX_start;     % update the total number of samples

LPF_ph = designfilt('lowpassfir', ...
    'PassbandFrequency',0.02,'StopbandFrequency',0.04, ...
    'PassbandRipple',0.1,'StopbandAttenuation',20, ...
    'DesignMethod','equiripple');

N_frm_mov =  T_mov/Ts;
Beat_stg = zeros(350,N_frm_mov,length(scan_use));

for ii = 1:length(scan_use)
    scan_mode_use = scan_mode;
    if strcmp(scan_mode,'both')
        if ii<=4
            scan_mode_use = 'rect';
        else 
            scan_mode_use = 'cross';
        end
    end
    k_stg = scan_use(ii);
    % Read USRP RX file
    file_ref = sprintf('%s_stg%d',scan_mode_use,k_stg); % refrence of experiment file name 
    file_USRP = [file_addr,'/',file_ref];
    
    file_name = [file_USRP,'_0.dat']; % I channel file name
    Beat_I = read_complex_binary2(file_name,N_spl,RX_start);      % read I channel file
    file_name = [file_USRP,'_1.dat'];     % Q channel file name 
    Beat_Q = read_complex_binary2(file_name,N_spl,RX_start);  % read Q channel file

    % IQ
    Beat_C = real(Beat_I) - 1j*real(Beat_Q);    % combine I Q channels to get complex signal 
    clear Beat_I Beat_Q;

    file_name = [file_USRP,'_2.dat'];     % Ref channel file name 
    Beat_ref = read_complex_binary2(file_name,N_spl_frm,RX_start);  % read Q channel file
    Beat_ref = real(Beat_ref);
    
    % Align ramp
    align_ramp_idx = sync_ramp_v4(Beat_ref(1:N_spl_frm),0);
    
    N_frm = floor((N_spl-align_ramp_idx+1)/N_spl_frm);    % update number of frames
    N_spl = N_frm*N_spl_frm;
    Beat_C = Beat_C(align_ramp_idx:(align_ramp_idx+N_spl-1));     % update RX data
    clear Beat_ref;
    
    if IF.DC_blk
        figure();
        subplot(2,1,1);
        plot(real(Beat_C(1:N_spl_frm))); hold on;
        subplot(2,1,2);
        plot(imag(Beat_C(1:N_spl_frm))); hold on;
        dcblker = dsp.DCBlocker;
        Beat_C = dcblker(Beat_C);
        subplot(2,1,1);
        plot(real(Beat_C(1:N_spl_frm)));
        subplot(2,1,2);
        plot(imag(Beat_C(1:N_spl_frm)));
        title('DC block');
    end
    
    % Normalize IQ channels 
    if IF.normal
        Beat_real_norm = mean(norm(real(Beat_C)),2);
        Beat_imag_norm = mean(norm(imag(Beat_C)),2);
        Beat_C = real(Beat_C)+1j*imag(Beat_C)*Beat_real_norm/Beat_imag_norm;
    end
    
    % Reshape Frame
    Beat_C = reshape(Beat_C,N_spl_frm,N_frm);
    
    % Separate fast and slow ramp samples 
    Beat_fr = Beat_C(1:N_spl_fr,:);
    clear Beat_C;
    
    % Cut off low high freq in ramp
    N_lfreq_spl = 50;
    N_hfreq_spl = 100;
    Beat_fr = Beat_fr([(N_lfreq_spl+1):(N_spl_fr/2-N_hfreq_spl),(N_spl_fr/2+1+N_hfreq_spl):(N_spl_fr-N_lfreq_spl)],:);
    N_spl_fr_up = (N_spl_fr - (N_lfreq_spl + N_hfreq_spl)*2)/2;
    
    Beat_fr_up = Beat_fr(1:N_spl_fr_up,:);
    
    % depth FFT
    N_FFT = N_spl_fr_up;
    FFT_fr = fft(Beat_fr_up,N_FFT,1);
    f_axis_fr = linspace(0,(Rs-Rs/N_FFT),N_FFT);    d_axis_fr = f_axis_fr/As1*c/2;
    yax_idx = find((d_axis_fr>int_cbl_len)&(d_axis_fr<(dist_ulim+int_cbl_len)));
    d_axis = d_axis_fr(yax_idx)-int_cbl_len;
    
    FFT_fr_mag = abs(FFT_fr);
    [~,pk_idx] = max(mean(FFT_fr_mag,2));
   
    % select the moving frames according to FFT phase
    
    mov_start_idx = zeros(1,10);
    
    figure(); 
    for kyd = (pk_idx-5):(pk_idx+5)
        FFT_ph = phase(FFT_fr(kyd,:));
        Beat_ph = phase(Beat_fr_up(kyd,:));
        FFT_ph_lpf = filtfilt(LPF_ph,FFT_ph);
        Beat_ph_lpf = filtfilt(LPF_ph,Beat_ph);
%         FFT_ph_lpf_drv = diff(FFT_ph_lpf,2);
        if IF.debug 
%             plot(FFT_ph); hold on;   
            subplot(2,1,1);
            plot(FFT_ph_lpf); hold on;
            subplot(2,1,2);
            plot(Beat_ph_lpf); hold on;
    %         plot(FFT_ph_lpf_drv*100); hold on;
        end
        FFT_ph_lpf_avg = mean(FFT_ph_lpf(1:300),2);
        mov_start_idx(kyd) = find(abs(FFT_ph_lpf-FFT_ph_lpf_avg)>abs(FFT_ph_lpf_avg)*0.1,1)-10;
    end
    mov_start_idx = mov_start_idx(find(mov_start_idx>1));
    mov_start_idx = round(median(mov_start_idx));
%     mov_start_idx = 469;
    mov_end_idx = mov_start_idx + N_frm_mov -1;
    if IF.debug
        subplot(2,1,1);
        plot(repmat(mov_start_idx,1,21),-10:10,'r*'); hold on;
        plot(repmat(mov_end_idx,1,21),-10:10,'r*'); hold on;
        subplot(2,1,2);
        plot(repmat(mov_start_idx,1,21),-10:10,'r*'); hold on;
        plot(repmat(mov_end_idx,1,21),-10:10,'r*'); hold on;
    end
%     FFT_fr = FFT_fr(mov_start_idx:mov_end_idx,:);
    Beat_fr_up = Beat_fr_up(:,mov_start_idx:mov_end_idx);
    Beat_stg(:,:,ii) = Beat_fr_up;
end

clear Beat_fr Beat_fr_up FFT_fr FFT_ph FFT_ph_lpf;

%% AoA
theta=0:0.1:180;    N_ang = length(theta);
AoA_P = zeros(length(yax_idx), length(theta), length(scan_use));
% AoA_C = 

N_FFT = N_spl_fr_up;
FFT_stg = fft(Beat_stg,N_FFT,1);

if IF.load
    load(sprintf('~/Dropbox/Project/Imaging with ML/result/1015/o_%s_mc2.mat',obj_ref),'AoA_P');
else
    for ks = 1:length(scan_use)
        dl = v_ls*Ts;
        lambda = c/((61-1.5/2)*1e9);

        for ii=1:length(theta)
            Vec = exp(-1j*2*pi*dl*(1:N_frm_mov)*cos(theta(ii)*pi/180)/lambda);
            VecR = repmat(Vec,length(yax_idx),1);
            AoA_P(:,ii,ks)=abs(sum(FFT_stg(yax_idx,:,ks).*VecR,2)).^2;    
        end
        clear Vec VecR;
    end
    
    if IF.save
        save(sprintf('~/Dropbox/Project/Imaging with ML/result/1015/o_%s_mc2.mat',obj_ref),'AoA_P');
    end
end

for ks = 1:length(scan_use)
    figure();
    h = surf(theta,d_axis,AoA_P(:,:,ks));
    ylabel('distance [m]');    xlabel('theta');
    title(['scan ',num2str(ks)],'Interpreter', 'none');    colormap jet
    caxis([0,clim]);    colorbar;
    set(h,'LineStyle','none');  view(2);
%     [Cx,Cy,AoA_C] = P2C(d_axis, theta, AoA_P(:,:,ks).', 0.04, clim);
end

%% Multiplying horizontal with vertical

H_P = repmat(AoA_P(:,:,2),[1,1,N_ang]);
V_P = repmat(AoA_P(:,:,1),[1,1,N_ang]);
V_P = permute(V_P,[1 3 2]);

sp_P = H_P.*V_P;

Nr = 70;
Na = 1801;

theta_ra = theta/180*pi;
theta_re = pi/2-theta/180*pi;

%% Background subtraction
if IF.save_bg
    sp_P_bg = sp_P;
    save(['~/Jayden/Imaging_ML/Matlab/matfile/multi_bg_',date_ref],'sp_P_bg');
elseif IF.bg_sub
    load(['~/Jayden/Imaging_ML/Matlab/matfile/multi_bg_',date_ref]);
    sp_P_sub = sp_P - sp_P_bg;
end

%% Projection
% filter by r
sp_P_near = sp_P(find((d_axis>2)&(d_axis<3)),:,:);
sp_P_near_proj = squeeze(max(sp_P_near,[],1));
figure();
h = surf(theta_re,theta_ra,sp_P_near_proj);
ylabel('theta');    xlabel('phi');
title('Current Near Projection 2-3m','Interpreter', 'none');    colormap jet
caxis([0,0.02e15]);  
colorbar;
set(h,'LineStyle','none');  view([90 90]);

sp_P_near = sp_P_bg(find((d_axis>2)&(d_axis<3)),:,:);
sp_P_near_proj = squeeze(max(sp_P_near,[],1));
figure();
h = surf(theta_re,theta_ra,sp_P_near_proj);
ylabel('theta');    xlabel('phi');
title('Background Near Projection 2-3m','Interpreter', 'none');    colormap jet
caxis([0,0.02e15]);  
colorbar;
set(h,'LineStyle','none');  view(2);

sp_P_near = sp_P_sub(find((d_axis>2)&(d_axis<3)),:,:);
sp_P_near_proj = squeeze(max(sp_P_near,[],1));
figure();
h = surf(theta_re,theta_ra,sp_P_near_proj);
ylabel('theta');    xlabel('phi');
title('BG Subtraction Near Projection 2-3m','Interpreter', 'none');    colormap jet
caxis([0,0.02e15]);  
colorbar;
set(h,'LineStyle','none');  view(2);

%% 3D Scatter
sp_P_th_u = 0.1e16;%max(sp_P(:))*0.01;
sp_P_ls = find(sp_P>(sp_P_th_u));
sp_P_N = length(sp_P_ls);

p_r =  mod(sp_P_ls,Nr);
p_r(find(p_r == 0)) = Nr; 

p_az =  mod(ceil(sp_P_ls/Nr),Na);
p_az(find(p_az == 0)) = Na;

p_el = ceil(sp_P_ls/(Nr*Na));

[px,py,pz] = sph2cart(theta_ra(p_az),theta_re(p_el),d_axis(p_r));
figure();
scatter3(px,py,pz); hold on;
xlabel('x');
ylabel('y');
zlabel('z');

% 
% if IF.bg_sub
%     sp_P_th_u = 0.01e16;%max(sp_P(:))*0.01;
%     sp_P_ls = find(sp_P_bg>(sp_P_th_u));
%     sp_P_N = length(sp_P_ls);
% 
%     p_r =  mod(sp_P_ls,Nr);
%     p_r(find(p_r == 0)) = Nr; 
% 
%     p_az =  mod(ceil(sp_P_ls/Nr),Na);
%     p_az(find(p_az == 0)) = Na;
% 
%     p_el = ceil(sp_P_ls/(Nr*Na));
% 
%     theta_ra = theta/180*pi;
%     theta_re = pi/2-theta/180*pi;
%     [px,py,pz] = sph2cart(theta_ra(p_az),theta_re(p_el),d_axis(p_r));
%     figure();
%     scatter3(px,py,pz,'r'); hold on;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     
%     sp_P_th_u = 0.01e16;%max(sp_P(:))*0.01;
%     sp_P_ls = find(sp_P_sub>(sp_P_th_u));
%     sp_P_N = length(sp_P_ls);
% 
%     p_r =  mod(sp_P_ls,Nr);
%     p_r(find(p_r == 0)) = Nr; 
% 
%     p_az =  mod(ceil(sp_P_ls/Nr),Na);
%     p_az(find(p_az == 0)) = Na;
% 
%     p_el = ceil(sp_P_ls/(Nr*Na));
% 
%     theta_ra = theta/180*pi;
%     theta_re = pi/2-theta/180*pi;
%     [px,py,pz] = sph2cart(theta_ra(p_az),theta_re(p_el),d_axis(p_r));
%     figure();
%     scatter3(px,py,pz,'r'); hold on;
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
% end

%%

% sp_P_th_d = sp_P_th_u/2;%max(sp_P(:))*0.01;
% sp_P_ls = find((sp_P<(sp_P_th_u))&(sp_P>(sp_P_th_d)));
% sp_P_N = length(sp_P_ls);
% 
% p_r =  mod(sp_P_ls,Nr);
% p_r(find(p_r == 0)) = Nr; 
% 
% p_az =  mod(ceil(sp_P_ls/Nr),Na);
% p_az(find(p_az == 0)) = Na;
% 
% p_el = ceil(sp_P_ls/(Nr*Na));
% 
% theta_ra = theta/180*pi;
% theta_re = pi/2-theta/180*pi;
% [px,py,pz] = sph2cart(theta_ra(p_az),theta_re(p_el),d_axis(p_r));
% scatter3(px,py,pz,'b');
% xlabel('x');
% ylabel('y');
% zlabel('z');


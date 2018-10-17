close all;
clc; clear;
addpath('../function');

%% Processing selection
IF = struct('debug',0, ...
            'align_ramp',1, ...
            'DC_blk',0, ...
            'normal',1, ...
            'bg_sub',1, ...
            'slow_ramp',0, ...
            'cut_lh_freq',1, ...
            'save_bg',0, ...
            'plt_phase',0, ...
            'cart',1, ...
            'save_plt',0);

%% Parameters
wvform = 'fsr_tri_6ms';

Rs = 5e5;   % USRP sampling rate 
T_mov = 30;
T_meas = 35;            % time of sampling in seconds
v_ls = 1.25e-2;
N_spl = T_meas * Rs;    % number of samples
inter_dist = 2.76;

date_ref = '1009';      % reference of experiment date
obj_ref = '1';
run_ref = 1;
% scan_mode = 'cross';
% k_stg = 4;
scan_mode = 'rect';
k_stg = 1;
ant_ref = 'o'; % 'o';
variable_library;
dist_ulim = 10;

clim = 100000000;

variable_library;
dist_ulim = 6;

%% Read USRP RX file
%     file_addr = sprintf('~/Dropbox/Project/Imaging with ML/data/%s/o%s_r%d',date_ref,obj_ref,run_ref);
file_addr = sprintf('~/Jayden/Imaging_ML/data/%s/o%s_r%d',date_ref,obj_ref,run_ref);

file_ref = sprintf('%s_stg%d',scan_mode,k_stg); % refrence of experiment file name 
file_USRP = [file_addr,'/',file_ref];

% Cut off USRP setup
RX_start = Rs/10;    % skip the first half second samples due to USRP start up
N_spl = N_spl-RX_start;     % update the total number of samples

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


%% Pre-processing
%Align ramp
if IF.align_ramp
    align_ramp_idx = sync_ramp_v4(Beat_ref(1:N_spl_frm),0);
else
    load('./matfile/RX_ref_05M.mat');
    sync_ramp_ar = zeros(N_spl_fr+N_spl_sr,1);
    for kc = 1:(N_spl_fr+N_spl_sr)
        sync_ramp_corr = corrcoef(RX_ref(kc:(kc+frm_len-1)),RX_ref_05M);
        sync_ramp_ar(kc)= sync_ramp_corr(2,1);
    end
    [~,align_ramp_idx] = max(sync_ramp_ar);
end

if IF.debug
    figure();
    plot(Beat_ref(1:N_spl_frm)); hold on;
    plot(align_ramp_idx,Beat_ref(align_ramp_idx), '*r');
    title(['Align ref ',file_ref],'Interpreter', 'none');
end

N_frm = floor((N_spl-align_ramp_idx+1)/N_spl_frm);    % update number of frames
N_spl = N_frm*N_spl_frm;
Beat_C = Beat_C(align_ramp_idx:(align_ramp_idx+N_spl-1));     % update RX data
clear Beat_ref;

if IF.debug
    figure();
    subplot(2,2,1); plot(real(Beat_C(1:N_spl_frm)));
    subplot(2,2,2); plot(imag(Beat_C(1:N_spl_frm)));
    subplot(2,2,3); plot(abs(Beat_C(1:N_spl_frm)));
    subplot(2,2,4); plot(angle(Beat_C(1:N_spl_frm)));
    sgtitle('Compare IQ channel spl')
end

% Dc removal filter
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

% Data loss
if IF.debug
    figure();
    plot(real(Beat_C(:)));
    title('received signal for loss detection');
end

% Reshape Frame
Beat_C = reshape(Beat_C,N_spl_frm,N_frm);

%% Separate fast and slow ramp samples 
Beat_fr = Beat_C(1:N_spl_fr,:);
if IF.slow_ramp
    Beat_sr = Beat_C(N_spl_fr+1:N_spl_frm,:);
end
clear Beat_C;

%% Cut off low high freq in ramp
N_lfreq_spl = 10;
N_hfreq_spl = 10;
Beat_fr = Beat_fr([(N_lfreq_spl+1):(N_spl_fr/2-N_hfreq_spl),(N_spl_fr/2+1+N_hfreq_spl):(N_spl_fr-N_lfreq_spl)],:);
N_spl_fr = N_spl_fr - (N_lfreq_spl + N_hfreq_spl)*2;
if IF.slow_ramp
    N_lfreq_spl = N_lfreq_spl*2;   N_hfreq_spl = N_hfreq_spl*2;
    Beat_sr =  Beat_fr([(N_lfreq_spl+1):(N_spl_fr/2-N_hfreq_spl),(N_spl_fr/2+1+N_hfreq_spl):(N_spl_fr-N_lfreq_spl)],:);
    N_spl_sr = N_spl_sr - (N_lfreq_spl + N_hfreq_spl)*2;
end

%% depth FFT
N_FFT = N_spl_fr;
FFT_fr = fft(Beat_fr,N_FFT,1);
f_axis_fr = linspace(0,(Rs-Rs/N_FFT),N_FFT);    d_axis_fr = f_axis_fr/As1*c/2;
yax_idx = find((d_axis_fr>inter_dist)&(d_axis_fr<(dist_ulim+inter_dist)));
d_axis = d_axis_fr(yax_idx)-inter_dist;

if IF.slow_ramp
    N_FFT = N_spl_sr;
    FFT_sr = fft(Beat_sr,N_FFT,1);
    f_axis_sr = linspace(0,(Rs-Rs/N_FFT),N_FFT);
    d_axis_sr = f_axis_fr/As2*c/2;
    d_axis = d_axis_sr(yax_idx);
end

%% depth FFT color map 
FFT_mag_fr = abs(FFT_fr);
figure();
h=surf(1:N_frm,d_axis,FFT_mag_fr(yax_idx,:));
title('Time domain with direct path cutoff');
xlabel('frame'); ylabel('distance [m]');
colormap jet;   caxis([0,30]);  colorbar;
set(h,'LineStyle','none');  view(2);

figure();
h=surf(1:N_frm,d_axis_fr,FFT_mag_fr);
title('Time domain start at 0');
xlabel('frame'); ylabel('distance [m]');
colormap jet;   caxis([0,30]);  colorbar;
set(h,'LineStyle','none');  view(2);
% clear FFT_mag_fr;

%% select the moving frames according to FFT phase
lin_l = 1000;
switch ant_ref 
    case'3'
        lin_r = 1500;
    case'12'
        lin_r = 1500;
    case 'o' 
        lin_r = 2000;
end

FFT_ph_dir = phase(FFT_fr(yax_idx(1),:));
% if IF.debug
    figure(); plot(FFT_ph_dir); hold on;
% end
FFT_ph_dir_slope = (FFT_ph_dir(lin_r)-FFT_ph_dir(lin_l))/(lin_r-lin_l);
mov_start_idx = round((lin_l+lin_r)/2-median(FFT_ph_dir((lin_l+1):lin_r))/FFT_ph_dir_slope);
N_frm_mov =  T_mov/Ts;

mov_start_idx = 558;

mov_end_idx = mov_start_idx + N_frm_mov -1;
clear Beat_ph_1
% if IF.debug
    plot(mov_start_idx,FFT_ph_dir(mov_start_idx),'r*'); 
    plot(mov_end_idx,FFT_ph_dir(mov_end_idx),'r*');
% end
FFT_fr = FFT_fr(:, mov_start_idx:mov_end_idx);

%% select the moving frames according to angle
% Beat_ph_1 = mean(phase(Beat_C),1);
% figure(); plot(Beat_ph_1); hold on;
% filt_lp = designfilt('lowpassfir', 'PassbandFrequency',0.05, ...
%     'StopbandFrequency',0.1,'DesignMethod','equiripple');
% Beat_ph_1 = filtfilt(filt_lp,Beat_ph_1);
% plot(Beat_ph_1);
% [~, mov_start_idx] = min(Beat_ph_1); mov_start_idx = mov_start_idx + 11;
% N_frm_mov =  T_mov/Ts;
% mov_end_idx = mov_start_idx + N_frm_mov -1;
% clear Beat_ph_1
% 
% Beat_C = Beat_C(:, mov_start_idx:mov_end_idx);

%% AoA
theta=0:0.1:180;
dl = v_ls*Ts;
lambda = c/((61-1.5/2)*1e9);

AoA_P = zeros(length(yax_idx), length(theta));
for ii=1:length(theta)
    Vec = exp(-1j*2*pi*dl*(1:N_frm_mov)*cos(theta(ii)*pi/180)/lambda);
    VecR = repmat(Vec,length(yax_idx),1);
    AoA_P(:,ii)=abs(sum(FFT_fr(yax_idx,:).*VecR,2)).^2;    
end

%theta = flip(theta);

figure();
h=surf(theta,d_axis,AoA_P);
% imagesc(theta,d_axis_fr(yax_idx),P(yax_idx,:));
xlabel('theta');    ylabel('distance [m]');
title(['current ',file_ref],'Interpreter', 'none');    colormap jet
caxis([0,clim]);    colorbar;
set(h,'LineStyle','none');  view(2);

if IF.save_bg
    AoA_bg = AoA_P;
    save(['../matfile/AoA_bg_',date_ref],'AoA_bg');
end

%% Plot
theta_plt = flip(theta);
figure('units','normalized','outerposition',[0 0 1 1]);
fig_polar_cu = polarPcolor(d_axis,theta_plt,AoA_P);
camroll(90);
title(['current polar ',file_ref],'Interpreter', 'none');    colormap jet
caxis([0,clim]); 

%% background subtraction
if IF.bg_sub
    load(['../matfile/AoA_bg_',date_ref]);
    % do nothing see what happens
    AoA_bg_sub = AoA_P - AoA_bg;
    
%     figure();
%     h=surf(theta,d_axis,P_bg);
%     % imagesc(theta,d_axis_fr(yax_idx),P(yax_idx,:));
%     xlabel('theta');    ylabel('distance [m]');
%     title('background');    colormap jet
%     caxis([0,50000000]);    colorbar;
%     set(h,'LineStyle','none');  view(2);
% 
%     figure();
%     h=surf(theta,d_axis,P_bg_sub);
%     % imagesc(theta,d_axis_fr(yax_idx),P(yax_idx,:));
%     xlabel('theta');    ylabel('distance [m]');
%     title(['background subtraction ',file_ref],'Interpreter', 'none');    colormap jet
%     caxis([-25000000,25000000]);    colorbar;
%     set(h,'LineStyle','none');  view(2);

    figure();
    polarPcolor(d_axis,theta_plt,AoA_bg);
    camroll(90);
    title(['back ground polar ',file_ref],'Interpreter', 'none');    colormap jet
    caxis([0,clim]); 
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    fig_polar_bg = polarPcolor(d_axis,theta_plt,AoA_bg_sub);
    camroll(90);
    title(['background subtraction ',file_ref],'Interpreter', 'none');    colormap jet
    caxis([-clim/2,clim/2]);
end

%% Cartesian
if IF.cart
    figure();
    [Cx,Cy,AoA_C] = P2C(d_axis, theta, AoA_P, 0.04, clim);
    [Cx, Cy] = meshgrid(Cx, Cy);
    fig_cart_cu = surf(Cx, Cy, AoA_C');
    title(['current cartitian ',file_ref],'Interpreter', 'none');    colormap jet
    caxis([0,clim]);
    view(2);
end

%% Save plots
if IF.save_plt
    saveas(fig_polar_cu,sprintf('~/Jayden/Imaging_ML/figure/%s/o%s_r%d_pcu.fig',date_ref,obj_ref,run_ref));
    saveas(fig_polar_cu,sprintf('~/Jayden/Imaging_ML/figure/%s/o%s_r%d_pcu.jpg',date_ref,obj_ref,run_ref));
    if IF.bg_sub
        saveas(fig_polar_bg,sprintf('~/Jayden/Imaging_ML/figure/%s/o%s_r%d_pbg.fig',date_ref,obj_ref,run_ref));
        saveas(fig_polar_bg,sprintf('~/Jayden/Imaging_ML/figure/%s/o%s_r%d_pbg.jpg',date_ref,obj_ref,run_ref));
    end
    if IF.cart
        saveas(fig_cart_cu,sprintf('~/Jayden/Imaging_ML/figure/%s/o%s_r%d_ccu.fig',date_ref,obj_ref,run_ref));
        saveas(fig_cart_cu,sprintf('~/Jayden/Imaging_ML/figure/%s/o%s_r%d_ccu.jpg',date_ref,obj_ref,run_ref));
    end
end

%% Phase
if IF.plt_phase
    figure();
    subplot(2,1,1);
    plot(phase(Beat_fr(1,:)));    hold on;
    subplot(2,1,2);
    plot(phase(Beat_sr(1,:)));    hold on;
    sgtitle('phase of the 1st sample');

    figure();
    subplot(2,1,1);
    plot(phase(Beat_fr(251,:)));    hold on;
    subplot(2,1,2);
    plot(phase(Beat_sr(501,:)));    hold on;
    sgtitle('phase of the middle sample');
end


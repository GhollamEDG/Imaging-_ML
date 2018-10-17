% multi_AoA_3D
close all;
clc;    clear all;
addpath('../function');

%% Processing selection
IF = struct('debug',1, ...
            'DC_blk',0, ...
            'normal',1, ...
            'bg_sub',1, ...
            'slow_ramp',0, ...
            'cut_lh_freq',1, ...
            'save',1,...
            'load',0);

%% Parameters
scan_mode = 'cross';%'rect';
wvform = 'fsr_tri_6ms';
ant_ref = 'd12';

date_ref = '1009';      % reference of experiment date
obj_ref = '1';
run_ref = 1;

Rs = 5e5;   % USRP sampling rate theta
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
%     
%     figure();
%     imagesc(d_axis_fr,1:length(FFT_fr),abs(FFT_fr));
    
%     FFT_fr_mag = abs(FFT_fr(:,1:find(d_axis_fr>int_cbl_len+1,1)));
%     [~,dir_idx] = max(mean(FFT_fr_mag,1));
   
    % select the moving frames according to FFT phase
    
    mov_start_idx = zeros(1,10);
    
    figure(); 
    for kyd = 1:round(length(yax_idx)/2)
        FFT_ph = phase(FFT_fr(yax_idx(kyd),:));
        FFT_ph_lpf = filtfilt(LPF_ph,FFT_ph);
%         FFT_ph_lpf_drv = diff(FFT_ph_lpf,2);
        if IF.debug 
%             plot(FFT_ph); hold on;
            plot(FFT_ph_lpf); hold on;
    %         plot(FFT_ph_lpf_drv*100); hold on;
        end
        FFT_ph_lpf_avg = mean(FFT_ph_lpf(1:300),2);
        mov_start_idx(kyd) = find(abs(FFT_ph_lpf-FFT_ph_lpf_avg)>abs(FFT_ph_lpf_avg)*0.2,1);
    end
    mov_start_idx = mov_start_idx(find(mov_start_idx>1));
    mov_start_idx = round(median(mov_start_idx));
    mov_start_idx = 386;
    mov_end_idx = mov_start_idx + N_frm_mov -1;
    if IF.debug
        plot(repmat(mov_start_idx,1,21),-10:10,'r*'); hold on;
        plot(repmat(mov_end_idx,1,21),-10:10,'r*'); hold on;
    end
%     FFT_fr = FFT_fr(mov_start_idx:mov_end_idx,:);
    Beat_fr_up = Beat_fr_up(:,mov_start_idx:mov_end_idx);
    Beat_stg(:,:,ii) = Beat_fr_up;
end

clear Beat_fr Beat_fr_up FFT_fr FFT_ph FFT_ph_lpf;

%% AoA
theta = (0:1:180)/180*pi;    N_az = length(theta);
phi = (45:1:135)/180*pi;     N_el = length(phi);

sp_P = zeros(length(yax_idx), N_az*N_el);

N_FFT = N_spl_fr_up; N_dist = length(yax_idx);
FFT_stg = fft(Beat_stg,N_FFT,1);
FFT_stg = reshape(FFT_stg(yax_idx,:,:),N_dist,[]);

dl = v_ls*Ts;
lambda = c/((61-1.5/2)*1e9);
array_idx = -N_frm_mov/2+1:N_frm_mov/2;
    
for ka = 1:N_el*N_az
    kaz =  mod(ka,N_az);
    if (kaz == 0)
        kaz = N_az;
    end
    kel = ceil(ka/N_az);
    
    cos_theta = cos(theta(kaz));
    sin_theta = sin(theta(kaz));
    cos_phi = cos(phi(kel));
    sin_phi = sin(phi(kel));
    cost_sinp = cos_theta * sin_phi;
    
%     Ph1 = 2*pi*dl*(array_idx*cost_sinp-N_frm_mov/2*cos_phi)/lambda;
%     Ph2 = 2*pi*dl*(array_idx*cos_phi+N_frm_mov/2*cost_sinp)/lambda;
%     Ph3 = 2*pi*dl*(flip(array_idx)*cost_sinp+N_frm_mov/2*cos_phi)/lambda;
%     Ph4 = 2*pi*dl*(flip(array_idx)*cos_phi-N_frm_mov/2*cost_sinp)/lambda;
%     Ph5 = 2*pi*dl*(array_idx*cos_phi)/lambda;
%     Ph6 = 2*pi*dl*(array_idx*cost_sinp)/lambda;

%     Ph5 = zeros(1,N_frm_mov);
%     Ph6 = Ph5;

%     Ph = [Ph1,Ph2,Ph3,Ph4,Ph5,Ph6];
    Ph = [Ph5,Ph6];
    Vec = exp(-1j*(Ph));
    VecR = repmat(Vec,N_dist,1);
    AoA_P = sum(FFT_stg.*VecR,2);
    sp_P(:,ka)=abs(AoA_P).^2;
end

if IF.save
    save(sprintf('~/Jayden/Imaging_ML/result/1015/o_%s_fa.mat',obj_ref),'sp_P');
elseif IF.load
    load(sprintf('~/Jayden/Imaging_ML/result/1015/o_%s_fa.mat',obj_ref),'sp_P');
end

sp_P = reshape(sp_P,N_dist,N_az,N_el);

%%
sp_P_th = max(sp_P(:))*0.1;
sp_P_ls = find(sp_P>(sp_P_th));
sp_P_N = length(sp_P_ls);

kr =  mod(sp_P_ls,N_dist);
kr(find(kr == 0)) = N_dist; 

kaz =  mod(ceil(sp_P_ls/N_dist),N_az);
kaz(find(kaz == 0)) = N_az;

kel = ceil(sp_P_ls/(N_dist*N_az));

phi_p = pi/2-phi;
[px,py,pz] = sph2cart(theta(kaz),phi_p(kel),d_axis(kr));
figure();
scatter3(px,py,pz);
xlabel('x');
ylabel('y');
zlabel('z');
% space_P = zeros();


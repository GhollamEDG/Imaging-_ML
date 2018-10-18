close all;
clc;    clear;
addpath('../function');

%%
date_ref = '1009'; 
obj_ref = '26';
run_ref = 1;

IF_main = struct('pre',0, ...
    'sv_int',0, ...
    'AoA',0, ...
    'save_bg',1, ...
    'bg_sub',1);

data_addr = sprintf('~/Jayden/Imaging_ML/data/%s/o%s_r%d',date_ref,obj_ref,run_ref);
result_addr = '~/Jayden/Imaging_ML/result/mlt_crs';

%
var_lib_mltcrs;

%% Preprocessing
pre_file = sprintf('%s/inter/%s_o%s_r%d_beat.mat',result_addr,date_ref,obj_ref,run_ref);
if IF_main.pre
    mlt_crs_pre;
    if IF_main.sv_int
        save(pre_file,'Beat_stg','');
    end
else
    load(pre_file,'Beat_stg');
end

%% FFT
N_FFT = N_mfreq_spl;
FFT_stg = fft(Beat_stg,N_FFT,1);
FFT_stg = FFT_stg(yax_idx,:,:);

%% AoA
AoA_file = sprintf('%s/inter/%s_o%s_r%d_AoA.mat',result_addr,date_ref,obj_ref,run_ref);

if IF_main.AoA
    mlt_crs_AoA;
    if IF_main.sv_int
        save(AoA_file,'AoA_P','');
    end
else
    load(AoA_file,'AoA_P');
end

%% Multiplying horizontal with vertical
H_P = repmat(AoA_P(:,:,2),[1,1,N_agl]);
V_P = repmat(AoA_P(:,:,1),[1,1,N_agl]);
V_P = permute(V_P,[1 3 2]);

sp_P = H_P.*V_P;

theta_ra = theta/180*pi;
theta_re = pi/2-theta/180*pi;

%% Background subtraction
if IF_main.save_bg
    sp_P_bg = sp_P;
    save(sprintf('%s/inter/mlt_bg%s_%s',result_addr,date_ref,ant_ref),'sp_P_bg');
end
if IF_main.bg_sub
    load(sprintf('%s/inter/mlt_bg%s_%s',result_addr,date_ref,ant_ref),'sp_P_bg');
    sp_P_sub = sp_P - sp_P_bg;
end

%%
sph_prj_2D;

%%
plt_sph_scatter;


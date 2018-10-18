c = 3e8;    % speed of light
f_pu = 0;   % USRP passband frequency

switch wvform
    case 'fsr_tri_6ms'
        % 2 slopes, 1.5 GHz in 1ms and 2ms, triangular ramp, total Ts=6ms
        dt_stp = 0.025e-6;      % time step of ramps
        df_stp1 = 37503.242;    % freq step of ramp1
        N_stp1 = 40000;         % number of steps of ramp1
        N_spl_fr = N_stp1*dt_stp*Rs*2;
        df_stp2 = 18751.621;    % freq step of ramp2
        N_stp2 = 80000;         % number of steps of ramp2
        N_spl_sr = N_stp2*dt_stp*Rs*2;
        N_spl_frm = N_spl_fr + N_spl_sr;

        f_b1 = 0e6;             % ramp start freq
        BW_FMCW = df_stp1*(N_stp1-1);
        f_b2 = f_b1+BW_FMCW;    % ramp end freq
        Ts = dt_stp*(N_stp1+N_stp2)*2;    % time span of a FMCW ramp cycle
        As1 = (BW_FMCW)/(dt_stp*N_stp1);    % freq slope of ramp1
        As2 = (BW_FMCW)/(dt_stp*N_stp2);    % freq slope of ramp2
end

% load(sprintf('./matfile/exp_list_%s.mat',date_ref),'exp_para');

% N_z_stp = exp_para(exp_idx).z_rng/exp_para(exp_idx).z_stp;

% obj_ref = exp_para(exp_idx).obj;
% run_ref = exp_para(exp_idx).run;

% var_lib_3D

%%
c = 3e8;    % speed of light
f_pu = 0;   % USRP passband frequency

% wvform parameters
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

% RX TX locations
TX_pos = struct('x',-0.24,'y',-0.035,'z',-0.155);

% linear stage parameters
x_ls_rng = 150;
z_ls_rng = 150;
x_ls_rng_h = floor(x_ls_rng/2);
z_ls_rng_h = floor(z_ls_rng/2);

% data collection (array scanning) parameters
% switch scan_mode
%     case 'rect'
%         N_scan_stg = 4;
%         scan_stg = struct('Num',4,'x_stp',[x_ls_rng,0,-x_ls_rng,0],'z_stp',[0,z_ls_rng,0,-z_ls_rng]);
%         scan_use = [1,2,3,4];
%     case 'cross'
%         N_scan_stg = 6;
%         scan_stg = struct('Num',6,'x_stp',[0,-x_ls_rng_h,0,x_ls_rng,0,-x_ls_rng_h],'z_stp',[z_ls_rng,0,-z_ls_rng_h,0,-z_ls_rng_h,0]);
%         scan_use = [1,4];
%     case 'both'
%         N_scan_stg = 10;
%         scan_stg = struct('Num',6,'x_stp',[0,-x_ls_rng_h,0,x_ls_rng,0,-x_ls_rng_h],'z_stp',[z_ls_rng,0,-z_ls_rng_h,0,-z_ls_rng_h,0]);
%         scan_use = [1,2,3,4,1,4];
% end




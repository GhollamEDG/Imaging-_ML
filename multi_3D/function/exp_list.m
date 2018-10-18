%%
clear;

IF_load = 1;
IF_save = 1;

%%
date_ref = '1001';      % reference of experiment date

if (IF_load)
    load(sprintf('./matfile/exp_list_%s.mat',date_ref),'exp_para'); 
else
    % initiate exp list
    field1 = 'obj'; value1 = '0';
    field2 = 'run'; value2 = 0;
    field3 = 'x_rng';  value3 = 0;
    field4 = 'z_rng';  value4 = 0;
    field5 = 'x_stp';  value5 = 0;
    field6 = 'z_stp';  value6 = 0;
    exp_para = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);
end

%% file
obj_ref = '1';
run_ref = 16;

%% Linear stage
x_ls_range = 180;
z_ls_range = 0;
x_ls_opt = 'rotate'; % conti or sng
x_ls_stp = 0;
N_x_stp = 1;
z_ls_stp = 0;
if z_ls_stp ~= 0
    N_z_stp = z_ls_range/ z_ls_stp;
else
    N_z_stp = 1;
end

%% new exp_para
if IF_load
    new_idx = size(exp_para,2)+1;
else 
    new_idx = 1;
end
exp_para(new_idx) .obj = obj_ref;
exp_para(new_idx) .run = run_ref;
exp_para(new_idx) .x_rng = x_ls_range;
exp_para(new_idx) .z_rng = z_ls_range;
exp_para(new_idx) .x_stp = x_ls_stp;
exp_para(new_idx) .z_stp = z_ls_stp;

%%
if (IF_save)
    save(sprintf('./matfile/exp_list_%s.mat',date_ref),'exp_para'); 
end
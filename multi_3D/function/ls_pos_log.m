function ls_pos = ls_pos_log(command,ls_x_new,ls_z_new,wv_len_half)
switch command
    case 'initiate'
        ls_pos = struct('x',ls_x_new,'z',ls_z_new,'wv_len',wv_len_half,'x_m',0,'z_m',0);
        save('~/Jayden/Imaging_ML/Matlab/matfile/ls_pos.mat','ls_pos'); 
    case 'reset_origin'
        load('~/Jayden/Imaging_ML/Matlab/matfile/ls_pos.mat','ls_pos'); 
        ls_pos.x = 0;
        ls_pos.z = 0;
        ls_pos.x_m = 0;
        ls_pos.z_m = 0;
        save('~/Jayden/Imaging_ML/Matlab/matfile/ls_pos.mat','ls_pos'); 
    case 'read'
        load('~/Jayden/Imaging_ML/Matlab/matfile/ls_pos.mat','ls_pos'); 
        fprintf(sprintf('ls pos: x=%d, z=%d\n',ls_pos.x,ls_pos.z));
    case 'update'
        load('~/Jayden/Imaging_ML/Matlab/matfile/ls_pos.mat','ls_pos');
        ls_pos.x = ls_x_new;
        ls_pos.z = ls_z_new;
        ls_pos.x_m = ls_x_new*wv_len_half;
        ls_pos.z_m = ls_z_new*wv_len_half;
        save('~/Jayden/Imaging_ML/Matlab/matfile/ls_pos.mat','ls_pos'); 
end
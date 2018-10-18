function ls_pos = rect2cross(ls_x,ls_z,x_ls_rng_h,dir,x_ls_port) 
    ls_pos = ls_pos_log('read',0,0,0);
    
        if dir>0
            if (ls_pos.x ~= 0)|(ls_pos.z ~= 0)
                error('Error. \nInvalid initial position')
            else
                fprintf('switch from rect to cross\n');
            end 
        else
            if (ls_pos.x ~= x_ls_rng_h)|(ls_pos.z ~= 0)
                error('Error. \nInvalid initial position')
            else
                fprintf('switch from cross to rect\n');
            end
        end
        N_x_stp = x_ls_rng_h*dir;
        x_arduino=sprintf('sudo echo -ne "%d,0\n" > /dev/ttyACM%d',N_x_stp,x_ls_port);
        system(x_arduino);
        ls_x_new = ls_pos.x + N_x_stp; 
        ls_pos = ls_pos_log('update',ls_x_new,ls_z,0.005);
end
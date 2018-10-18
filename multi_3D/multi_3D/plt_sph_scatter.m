%% 3D Scatter
sp_P_th_u = 0.1e16;%max(sp_P(:))*0.01;
sp_P_ls = find(sp_P>(sp_P_th_u));
sp_P_N = length(sp_P_ls);

p_r =  mod(sp_P_ls,N_r);
p_r(find(p_r == 0)) = N_r; 

p_az =  mod(ceil(sp_P_ls/N_r),N_agl);
p_az(find(p_az == 0)) = N_agl;

p_el = ceil(sp_P_ls/(N_r*N_agl));

[px,py,pz] = sph2cart(theta_ra(p_az),theta_re(p_el),d_axis(p_r));
figure();
scatter3(px,py,pz); hold on;
xlabel('x');
ylabel('y');
zlabel('z');

% W
% if IF.bg_sub
%     sp_P_th_u = 0.01e16;%max(sp_P(:))*0.01;
%     sp_P_ls = find(sp_P_bg>(sp_P_th_u));
%     sp_P_N = length(sp_P_ls);
% 
%     p_r =  mod(sp_P_ls,N_r);
%     p_r(find(p_r == 0)) = N_r; 
% 
%     p_az =  mod(ceil(sp_P_ls/N_r),N_agl);
%     p_az(find(p_az == 0)) = N_agl;
% 
%     p_el = ceil(sp_P_ls/(N_r*N_agl));
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
%     p_r =  mod(sp_P_ls,N_r);
%     p_r(find(p_r == 0)) = N_r; 
% 
%     p_az =  mod(ceil(sp_P_ls/N_r),N_agl);
%     p_az(find(p_az == 0)) = N_agl;
% 
%     p_el = ceil(sp_P_ls/(N_r*N_agl));
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
% p_r =  mod(sp_P_ls,N_r);
% p_r(find(p_r == 0)) = N_r; 
% 
% p_az =  mod(ceil(sp_P_ls/N_r),N_agl);
% p_az(find(p_az == 0)) = N_agl;
% 
% p_el = ceil(sp_P_ls/(N_r*N_agl));
% 
% theta_ra = theta/180*pi;
% theta_re = pi/2-theta/180*pi;
% [px,py,pz] = sph2cart(theta_ra(p_az),theta_re(p_el),d_axis(p_r));
% scatter3(px,py,pz,'b');
% xlabel('x');
% ylabel('y');
% zlabel('z');


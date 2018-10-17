%% AoA
theta=linspace(0,180,N_agl);
AoA_P = zeros(N_r, N_agl, length(stg_use));

for ks = 1:length(stg_use)
    dl = v_ls*Ts;
    lambda = c/((61-1.5/2)*1e9);
    for ii=1:N_agl
        Vec = exp(-1j*2*pi*dl*(1:N_frm_mov)*cos(theta(ii)*pi/180)/lambda);
        VecR = repmat(Vec,N_r,1);
        AoA_P(:,ii,ks)=abs(sum(FFT_stg(:,:,ks).*VecR,2)).^2;    
    end
    clear Vec VecR;
end

for ks = 1:length(stg_use)
    figure();
    h = surf(theta,d_axis,AoA_P(:,:,ks));
    ylabel('distance [m]');    xlabel('theta');
    title(['scan ',num2str(ks)],'Interpreter', 'none');    colormap jet
    caxis([0,clim]);    colorbar;
    set(h,'LineStyle','none');  view(2);
%     [Cx,Cy,AoA_C] = P2C(d_axis, theta, AoA_P(:,:,ks).', 0.04, clim);
end
%% Projection
% filter by r
sp_P_near = sp_P(find((d_axis>2)&(d_axis<3)),:,:);
sp_P_near_proj = squeeze(max(sp_P_near,[],1));
figure();
h = surf(theta_re,theta_ra,sp_P_near_proj);
ylabel('theta');    xlabel('phi');
title('Current Near Projection 2-3m','Interpreter', 'none');    colormap jet
caxis([0,0.02e15]);  
colorbar;
set(h,'LineStyle','none');  view([90 90]);

if IF_main.bg_sub
    sp_P_near = sp_P_bg(find((d_axis>2)&(d_axis<3)),:,:);
    sp_P_near_proj = squeeze(max(sp_P_near,[],1));
    figure();
    h = surf(theta_re,theta_ra,sp_P_near_proj);
    ylabel('theta');    xlabel('phi');
    title('Background Near Projection 2-3m','Interpreter', 'none');    colormap jet
    caxis([0,0.02e15]);  
    colorbar;
    set(h,'LineStyle','none');  view(2);

    sp_P_near = sp_P_sub(find((d_axis>2)&(d_axis<3)),:,:);
    sp_P_near_proj = squeeze(max(sp_P_near,[],1));
    figure();
    h = surf(theta_re,theta_ra,sp_P_near_proj);
    ylabel('theta');    xlabel('phi');
    title('BG Subtraction Near Projection 2-3m','Interpreter', 'none');    colormap jet
    caxis([0,0.02e15]);  
    colorbar;
    set(h,'LineStyle','none');  view(2);
end


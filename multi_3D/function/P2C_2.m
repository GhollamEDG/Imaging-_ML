function [X, Y, Q] = P2C_2(r, theta, P, gran, c_lim)

    rmax = max(r);
    P = P';
    
    [R, Theta] = meshgrid(r, theta);
    [X1, Y1] = pol2cart(Theta*pi/180, R);
    S = scatteredInterpolant(X1(:), Y1(:), P(:), 'linear', 'none');
    
    [X,Y] = meshgrid(-rmax:gran:rmax, 0:gran:rmax);

    Q = S(X, Y);
   
    
    figure();
    surf(X, Y, Q);
    caxis([0,c_lim]); 
    grid off;
    colormap jet;    colorbar;
    view(2);
end

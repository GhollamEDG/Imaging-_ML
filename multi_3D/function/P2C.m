function [X, Y, Q] = P2C(r, theta, P, gran, c_lim)


% sohrab: this code interpolates polar into cartesian by simply taking the weighted average of
% the four polar points surrounding particular cartesian point

R = max(r);

    X = -R:gran:R;
    Y = 0:gran:R;

    Nx = length(X);
    Ny = length(Y);

    Q = zeros(length(X) , length(Y));

    for xx = 1:Nx
        for yy = 1:Ny

            theta_xy = 180/pi * atan(Y(yy)/X(xx));
            if(theta_xy < 0)
                theta_xy = 180 + theta_xy;
            end

           r_xy = sqrt(X(xx)^2 + Y(yy)^2);

           idx_theta = find(theta >= theta_xy, 1);
           if(idx_theta == 1)
               theta1 = idx_theta;
           else
               theta1 = idx_theta-1;
           end

           if(idx_theta == length(theta))
               theta2 = length(theta);
           else
               theta2 = idx_theta;
           end


           idx_r = find(r >= r_xy, 1);

           if(idx_r == 1)
               r1 = 1;
           else
               r1 = idx_r-1;
           end

           if(r_xy > R)
               Q(xx, yy) = 0;
               continue;
           else
              r2 = idx_r;
           end
           P1 = P(r1,theta1);
           P2 = P(r2,theta1);
           P3 = P(r1,theta2);
           P4 = P(r2,theta2);
           r1 = r(r1); r2 = r(r2);
           theta1 = theta(theta1);
           theta2 = theta(theta2);

           weight1 = 1/sqrt( (X(xx) - r1*cos(theta1) )^2 + ( Y(yy) - r1*sin(theta1) )^2);
           weight2 = 1/sqrt( (X(xx) - r2*cos(theta1) )^2 + ( Y(yy) - r2*sin(theta1) )^2);
           weight3 = 1/sqrt( (X(xx) - r1*cos(theta2) )^2 + ( Y(yy) - r1*sin(theta2) )^2);
           weight4 = 1/sqrt( (X(xx) - r2*cos(theta2) )^2 + ( Y(yy) - r2*sin(theta2) )^2);

           Q(xx, yy) = (P1*weight1 + P2*weight2 + ...
               P3*weight3 + P4*weight4) / (weight1 + weight2 + weight3 + weight4);


        end
    end
    %%

    % [MX, MY] = meshgrid(X, Y);
    % figure();
    % surf(MX, MY, Q');
    % colormap jet
    % caxis([0,c_lim]);
    % colormap jet;    colorbar;
    % view(2);

function [pos_pvt,pdop] = triangulate_generic(sns_pos, az_meas, el_meas, rng_meas, init_sol, iter_idx, triang_type)

% triang_type = 1; angles only
% triang_type = 2; ranges only
% triang_type = 3; angles + ranges
% plott3(init_sol(1:3),'s');
if nargin == 4
    triang_type = 1; % default angles triangulation
end
if ~ismember(triang_type,[1 2 3])
    triang_type = 1;
end


dim = 3;
num_sns = length(az_meas(:,1));

pos_now = init_sol(1:dim);
err = inf;
pdop = 100;

count_tmp = 0;
count_max = inf; %10000;

if num_sns >= 2
    while err > 1e-3  && count_tmp<=count_max
        count_tmp = count_tmp + 1;

        delta_x = repmat(pos_now(1),num_sns,1)-sns_pos(1,:)';
        delta_y = repmat(pos_now(2),num_sns,1)-sns_pos(2,:)';
        delta_z = repmat(pos_now(3),num_sns,1)-sns_pos(3,:)';
        
        range_est = sqrt(delta_x.^2+delta_y.^2+delta_z.^2);
        range_g_est = sqrt(delta_x.^2+delta_y.^2);
        
        dh1dx = -delta_y ./ range_g_est.^2;
        dh1dy = delta_x ./ range_g_est.^2;
        dh1dz = zeros(num_sns,1);
        
        dh2dx = -delta_x .* delta_z ./ range_g_est ./ range_est.^2;
        dh2dy = -delta_y .* delta_z ./ range_g_est ./ range_est.^2;
        dh2dz = range_g_est ./ range_est.^2;
        
        dh3dx = delta_x ./ range_est;
        dh3dy = delta_y ./ range_est;
        dh3dz = delta_z ./ range_est;
        
        az_comp = atan2(delta_y,delta_x);
        az_comp(az_comp<0) = az_comp(az_comp<0)+2*pi;
        el_comp = atan2(delta_z,range_g_est);
        rng_comp = range_est;
        
        rng_err = rng_meas - rng_comp;
        el_err = el_meas - el_comp;
        az_err = az_meas - az_comp;
        
%         el_err = mod(el_err,2*pi);
%         el_err(el_err>pi) = el_err(el_err>pi)-2*pi;
%         az_err = mod(az_err,2*pi);
%         az_err(az_err>pi) = az_err(az_err>pi)-2*pi;
        % proper angle fixing
        if any(el_err > 3*pi/2)
            el_err(el_err > 3*pi/2) = el_err(el_err > 3*pi/2) - 2*pi;
        end
        if any(el_err < -3*pi/2)
            el_err(el_err < -3*pi/2) = el_err(el_err < -3*pi/2) + 2*pi;
        end
        if any(el_err > pi/2)
            el_err(el_err > pi/2) = el_err(el_err > pi/2) - pi;
        end
        if any(el_err < -pi/2)
            el_err(el_err < -pi/2) = el_err(el_err < -pi/2) + pi;
        end
        %--------------------
        if any(az_err > 3*pi/2)
            az_err(az_err > 3*pi/2) = az_err(az_err > 3*pi/2) - 2*pi;
        end
        if any(az_err < -3*pi/2)
            az_err(az_err < -3*pi/2) = az_err(az_err < -3*pi/2) + 2*pi;
        end
        if any(az_err > pi/2)
            az_err(az_err > pi/2) = az_err(az_err > pi/2) - pi;
        end
        if any(az_err < -pi/2)
            az_err(az_err < -pi/2) = az_err(az_err < -pi/2) + pi;
        end
        % proper angle fixing
        
                
        switch triang_type
            case 1 % azimuth, elevation
                ang_err = [az_err;el_err];
                Dpr_now = [dh1dx dh1dy dh1dz; dh2dx dh2dy dh2dz];
            case 2 % range
                ang_err = [rng_err];
                Dpr_now = [dh3dx dh3dy dh3dz];
            case 3 % azimuth, elevation, range
                ang_err = [az_err;el_err;rng_err];
                Dpr_now = [dh1dx dh1dy dh1dz; dh2dx dh2dy dh2dz; dh3dx dh3dy dh3dz];
        end
        
        new_err_vec = Dpr_now\ang_err;
        
        pos_now = pos_now + new_err_vec;
%         plott3(pos_now,'o')
        [pos_now'];
        if norm(new_err_vec)>1.1*err
            disp(['triagulation error is increasing ...  (', num2str(iter_idx), ')']);
        end
        err = norm(new_err_vec);
    end
    %     diag_Dpr = diag(inv(Dpr_now'*Dpr_now));
    diag_Dpr = diag(eye(3)/(Dpr_now'*Dpr_now));
    pdop = sqrt(sum(diag_Dpr(1:dim)));
else
    pos_pvt = init_sol;
end

pos_pvt = pos_now;

return




function [datas, pdops, SphMat, triang_sol] = TriangulateManyGeneric(az_mat, el_mat, rng_mat, sns_pos, triang_type)
% if triang_type==1
%     1. triangulate using (az,el) in order to calc (rng)
%     2. convert (az,el,rng) to (x,y,z), and put it into datas
%     3. GT (triang_sol) = average all (x,y,z) from sensors
%     
% if triang_type==3
%     1. convert (az,el,rng) to (x,y,z), and put it into datas
%     2. GT (triang_sol) = average all (x,y,z) from sensors

num_points = length(az_mat(:,1));
num_sns = length(sns_pos(1,:));
num_dim = 3;
track_triang = nan(num_dim,num_points);
pdops = nan(1,num_points);

rng_est = nan(num_points,num_sns);
az_est = nan(num_points,num_sns);
el_est = nan(num_points,num_sns);
datas = []; 


for ii=1:num_points
    
    az_meas = az_mat(ii,:)';
    el_meas = el_mat(ii,:)';
    rng_meas = rng_mat(ii,:)';
    if ii == 1
        init_sol = [0;0;0];
        if triang_type == 3
            [x_0,y_0,z_0] = sph2cart(az_mat(1,ii),el_mat(1,ii),rng_mat(1,ii));
            init_sol = [x_0;y_0;z_0]+sns_pos(:,1); % heuristics
        end
    else
        init_sol = track_triang(:,ii-1); % this is without any thinking just to make it run
    end
    if triang_type~=3
        [track_triang(:,ii),pdops(ii)] = triangulate_generic(sns_pos, az_meas, el_meas, rng_meas, init_sol, ii, triang_type);
    end
    %     plott3(track_triang(:,ii),'d','markerfacecolor','r')
    1;
end


for ii=1:num_sns
    
    if triang_type~=3
        
        delta_x = (track_triang(1,:)-repmat(sns_pos(1,ii),1,num_points))';
        delta_y = (track_triang(2,:)-repmat(sns_pos(2,ii),1,num_points))';
        delta_z = (track_triang(3,:)-repmat(sns_pos(3,ii),1,num_points))';
        
        
        %
        rng_comp = sqrt(delta_x.^2+delta_y.^2+delta_z.^2);
        rng_comp_ref = sqrt(sum((track_triang-repmat(sns_pos(:,ii),1,num_points)).^2))';
        range_g_now = sqrt(delta_x.^2+delta_y.^2);
        
        az_comp = atan2(delta_y,delta_x);
        az_comp(az_comp<0) = az_comp(az_comp<0)+2*pi;
        el_comp = atan2(delta_z,range_g_now);
        
        if triang_type==1 %(az,el) from sensor, we needed only to calc (rng)
            rng_est(:,ii) = rng_comp; %rng_comp_ref;
        end
        if triang_type==2 %(rng) from sensor, we needed only to calc (az,el)
            az_est(:,ii) = az_comp;
            el_est(:,ii) = el_comp;
        end
    end
    if triang_type == 1
        [x_t,y_t,z_t] = sph2cart(az_mat(:,ii),el_mat(:,ii),rng_est(:,ii));
        SphMat.az_mat(:,ii)  = az_mat(:,ii);
        SphMat.el_mat(:,ii)  = el_mat(:,ii);
        SphMat.rng_mat(:,ii) = rng_est(:,ii);
    elseif triang_type == 2
        [x_t,y_t,z_t] = sph2cart(az_est(:,ii),el_est(:,ii),rng_mat(:,ii));
        SphMat.az_mat(:,ii)  = az_est(:,ii);
        SphMat.el_mat(:,ii)  = el_est(:,ii);
        SphMat.rng_mat(:,ii) = rng_mat(:,ii);
    else
        [x_t,y_t,z_t] = sph2cart(az_mat(:,ii),el_mat(:,ii),rng_mat(:,ii));
        SphMat.az_mat(:,ii)  = az_mat(:,ii);
        SphMat.el_mat(:,ii)  = el_mat(:,ii);
        SphMat.rng_mat(:,ii) = rng_mat(:,ii);
    end
    %     [x_t,y_t,z_t] = sph2cart(az_mat(:,ii),el_mat(:,ii),rng_mat(:,ii));
    datas(ii).xyz = [x_t,y_t,z_t];
    
end

if (triang_type == 3) % track_triang is full of nans
    si=1;
    track_triang = datas(si).xyz + repmat(sns_pos(:,si),1,num_points)';
    for si=2:num_sns
        track_triang = track_triang + datas(si).xyz + repmat(sns_pos(:,si),1,num_points)';
    end
    track_triang = track_triang/num_sns;
end

triang_sol = track_triang;
if (triang_type ~= 3)
    triang_sol = triang_sol';
end
end
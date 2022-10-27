close all;
clear all;
clc
% randn('seed', 1);
% rand('seed', 1);

disp_flag = 1;

%% scenario
sensor_type = 1;1;3; %3; %
MC = 1;Noise = 0.2;
Samples = 90;
% for Samples = [3, 6, 9, 12, 15, 30, 90, 300, 900], %10, %
version_new = 1;
seg_num = 5;omega_sign = -1; ang_freq = 225; omega_deg = 1;
for Ns = 8;
    % Samples = 90; [30] %[90, 300, 900], %10, %
    init_states = [0,0,1000,100,0,10]';
    track = GenerateScenarioSynthetic_new(seg_num, init_states, omega_sign, ang_freq, omega_deg);
    time_interval = 1:length(track(1,:));
    
    NNNtracks = size(track,2);
    
    skip = floor(NNNtracks/(Samples-1)); %10;
    track = track(:,1:skip:end);
    times = time_interval(1:skip:end);
    num_points = length(track(1,:));
    
    %% sensors
    % circle
    sns_pos = ...
        [0     20000 30000  20000  0     -20000 -30000 -20000
        30000 20000 0     -20000 -30000 -20000  0      20000
        0     0     3000      0      0      0      0      0     ]/3*4;
       
    max_sens = numel(sns_pos(1,:));
    sens_inds = 1:Ns;%randperm(max_sens,Ns);
    sns_pos = sns_pos(:,sens_inds);
      
    sgm_az  = Noise*0.5*3e-3; %1.5:1.5:9[mRad]
    sgm_el  = Noise*0.5*3e-3; %1.5:1.5:9[mRad]
    sgm_rng = Noise*0.5*10;   %5:5:30[m]
    
    rot_errs_flag = 1;
   
    max_ang = 80;
    rot_psi   = (rand(1,Ns)*max_ang-max_ang/2);
    rot_theta = (rand(1,Ns)*max_ang-max_ang/2);
    rot_phi   = (rand(1,Ns)*max_ang-max_ang/2);
    
    rot_psi = rot_psi*pi/180;
    rot_theta = rot_theta*pi/180;
    rot_phi = rot_phi*pi/180;
    
    num_sns = size(sns_pos,2);
    
    MC_err = [];
    for kkkk=1:MC,
        
        %% GroundTruth
        GT.track     = track;
        GT.rot_psi   = rot_psi;
        GT.rot_theta = rot_theta;
        GT.rot_phi   = rot_phi;
        
        
        %% Calculate Sensor Measurements
        az_mat  = []; az_mat_nom  = [];
        el_mat  = []; el_mat_nom  = [];
        rng_mat = []; rng_mat_nom = [];
        for si = 1:num_sns
            tracks{si} = [];
            
            A_rot = EKF_Euler2DCM(rot_phi(si),rot_theta(si),rot_psi(si)); %angle2dcm(rot_phi(si),rot_theta(si),rot_psi(si));
            tracks{si}.track_rot = A_rot * (track - repmat(sns_pos(:,si),1,num_points)) + repmat(sns_pos(:,si),1,num_points);
            tracks{si}.track_nom = eye(3) * (track - repmat(sns_pos(:,si),1,num_points)) + repmat(sns_pos(:,si),1,num_points);
            
            delta_x_rot = tracks{si}.track_rot(1,:) - sns_pos(1,si);
            delta_y_rot = tracks{si}.track_rot(2,:) - sns_pos(2,si);
            delta_z_rot = tracks{si}.track_rot(3,:) - sns_pos(3,si);
            
            delta_x_nom = tracks{si}.track_nom(1,:) - sns_pos(1,si);
            delta_y_nom = tracks{si}.track_nom(2,:) - sns_pos(2,si);
            delta_z_nom = tracks{si}.track_nom(3,:) - sns_pos(3,si);
            
            tracks{si}.az_rot = atan2(delta_y_rot,delta_x_rot) + sgm_az * randn(size(delta_y_rot));
            tracks{si}.az_rot(tracks{si}.az_rot<0) = tracks{si}.az_rot(tracks{si}.az_rot<0)+2*pi;
            tracks{si}.el_rot = atan2(delta_z_rot, sqrt(delta_x_rot.^2+delta_y_rot.^2)) + sgm_el*randn(size(delta_z_rot));
            tracks{si}.range  = sqrt(delta_x_rot.^2 + delta_y_rot.^2 + delta_z_rot.^2) + sgm_rng*randn(size(delta_x_rot));
            
            tracks{si}.az_nom = atan2(delta_y_nom,delta_x_nom) + sgm_az * randn(size(delta_y_nom));
            tracks{si}.az_nom(tracks{si}.az_nom<0) = tracks{si}.az_nom(tracks{si}.az_nom<0)+2*pi;
            tracks{si}.el_nom = atan2(delta_z_nom, sqrt(delta_x_nom.^2+delta_y_nom.^2)) + sgm_el*randn(size(delta_z_nom));
            tracks{si}.range_nom  = sqrt(delta_x_nom.^2 + delta_y_nom.^2 + delta_z_nom.^2) + sgm_rng*randn(size(delta_x_nom));
            
            [x_t,y_t,z_t] = sph2cart(tracks{si}.az_rot, tracks{si}.el_rot, tracks{si}.range);
            x_noisy = x_t + sns_pos(1,si);
            y_noisy = y_t + sns_pos(2,si);
            z_noisy = z_t + sns_pos(3,si);
            tracks{si}.track_rot_noisy = [x_noisy;y_noisy;z_noisy];
            
            [x_tnom,y_tnom,z_tnom] = sph2cart(tracks{si}.az_nom, tracks{si}.el_nom, tracks{si}.range_nom);
            x_noisy_nom = x_tnom + sns_pos(1,si);
            y_noisy_nom = y_tnom + sns_pos(2,si);
            z_noisy_nom = z_tnom + sns_pos(3,si);
            tracks{si}.track_nom_noisy = [x_noisy_nom;y_noisy_nom;z_noisy_nom];
            
            
            az_mat  = [az_mat,  tracks{si}.az_rot'];
            el_mat  = [el_mat,  tracks{si}.el_rot'];
            rng_mat = [rng_mat, tracks{si}.range' ];
            
            az_mat_nom  = [az_mat_nom,  tracks{si}.az_nom'];
            el_mat_nom  = [el_mat_nom,  tracks{si}.el_nom'];
            rng_mat_nom = [rng_mat_nom, tracks{si}.range_nom'];
        end
        
        %--------------------------------------------------------------------------
        if disp_flag
            figure(100); hold on; grid on;
            % set(gcf,'Position',[-1800 528 560 420])
            %plott3(track, 'b', 'LineWidth', 2);
            plott3(track, 'b','linewidth',3);
            % plott3(track(:,1), 'sb', 'markerfacecolor', 'b');
            xlabel('X(m)');
            ylabel('Y(m)');
            zlabel('Z(m)');
            str_leg{1} = 'Ground Truth';
            for si = 1:num_sns
                plott3(tracks{si}.track_rot_noisy,'--','linewidth',1);
                str_leg{si+1} = ['Sensor ',num2str(si)];
            end
            % plott3(track(:,1),'bo');
            if MC>=1
                plott3(sns_pos,'mo','markerfacecolor','m');
                for ii=1:num_sns,
                    text(sns_pos(1,ii)+700,sns_pos(2,ii)+700,sns_pos(3,ii),num2str(ii),'Color',[1 0 1]);
                end
            end
            axis equal;
            legend(str_leg);
            %xlim([-22000 42000]);
        end
        
        [~, pdops] = TriangulateManyGeneric(az_mat_nom, el_mat_nom, rng_mat_nom, sns_pos, sensor_type);
        if disp_flag
            %             figure(300); hold on; grid on;
            %             plot(times,pdops);
            %             xlabel('Time(s)')
            %             ylabel('PDOP')
        end
        
        err_1 = []; err_2 = [];
        angles1 = []; angles2 = [];
        [err_1, angles1, datas1, err_GT] = CalibrateBias_universal(az_mat, el_mat, rng_mat, sns_pos, sensor_type, 0, GT, version_new);
        disp('GT_angles')
        disp(180/pi*[GT.rot_psi(1:Ns);GT.rot_theta(1:Ns);GT.rot_phi(1:Ns)])
        disp('Calc_angles')
        disp(angles1)
        disp('Err_angles [deg]')
        disp(abs(180/pi*[GT.rot_psi(1:Ns);GT.rot_theta(1:Ns);GT.rot_phi(1:Ns)]-angles1))
        disp('Err_angles [mRad]')
        disp(pi/0.180*abs(180/pi*[GT.rot_psi(1:Ns);GT.rot_theta(1:Ns);GT.rot_phi(1:Ns)]-angles1))

        
        if ~isempty(err_1) && ~isempty(err_2)
            if numel(err_1) > numel(err_2)
                err_2 = [err_2, err_2(end)*ones(1,numel(err_1)-numel(err_2))];
            else
                err_1 = [err_1, err_1(end)*ones(1,numel(err_2)-numel(err_1))];
            end
        end
        err_vec = [err_1;err_2]';
        
        if disp_flag
            GcfPos = get(gcf,'Position');
            GcfPos2 = GcfPos;
            % GcfPos2(1) = GcfPos2(1) + 600;
            
            figure(200); hold on; grid on
            set(gcf,'Position',GcfPos2)
            % plot(log10(err_vec), '.-');
            plot(err_vec, '.-');
            xlabel('Iteration number')
            ylabel('Normalized Error (m)')
        end
        
        % disp([[err_GT(:).psi]', [err_GT(:).theta]', [err_GT(:).phi]', [err_GT(:).Cart]']);
        disp([[err_GT(end).psi]', [err_GT(end).theta]', [err_GT(end).phi]', [err_GT(end).Cart]']);
        
        
        MC_err(kkkk,:) = [[err_GT(end).psi]', [err_GT(end).theta]', [err_GT(end).phi]', [err_GT(end).Cart]'];
        MC_Niter(kkkk) = numel(err_GT);
        disp(['MC #',num2str(kkkk)]);
        MC_err_GT{kkkk} = err_GT;
    end % MC
    
    disp(mean(MC_err))
    disp(std(MC_err))
%     disp(['MC2_sensors__N',num2str(num_sns),'_type',num2str(sensor_type),'_NumPoints',num2str(num_points),'.mat']);
%     save(['MC2_sensors__N',num2str(num_sns),'_type',num2str(sensor_type),'_NumPoints',num2str(num_points),'.mat'])
    
end

1;
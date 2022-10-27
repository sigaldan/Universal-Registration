function [err_vec, angles, datas, err_GT] = CalibrateBias_universal(az_mat, el_mat, rng_mat, sns_pos, sensor_type, plot_fig_100, GT, version_new)

% sensor_type = 1; angles only
% sensor_type = 2; ranges only
% sensor_type = 3; angles + ranges

% --------------
flag_sub_triang = 0; %0-once in a while, 1-after dual pair, 2-after each pair
flag_sub_plot = 0;
intermediate_plots = 1; %1; %
M = struct([]);
% --------------

num_sns   = length(sns_pos(1,:));
num_points = size(az_mat,1);

err = inf;
err_vec = [];
err_GT = [];

if sensor_type == 1
    is2D = 1;
else
    is2D = 0;
end

%% GroundTruth
GT_track     = GT.track;
GT_rot_psi   = GT.rot_psi;
GT_rot_theta = GT.rot_theta;
GT_rot_phi   = GT.rot_phi;


%%
[datas, dops, SphMat, track_triang] = TriangulateManyGeneric(az_mat, el_mat, rng_mat, sns_pos, sensor_type);
datas_original = datas;

%% disp initial errors
errs = [];
for si1 = 1:num_sns
    for si2 = si1+1:num_sns
        dif12_ORG2 = datas(si1).xyz' + repmat(sns_pos(:,si1)-sns_pos(:,si2),1,num_points) - datas(si2).xyz';
        erri = mean(sqrt(sum(dif12_ORG2.^2)));
        errs = [errs, erri];
    end
end
disp([' 0: ', num2str([errs, sum(errs)/((num_sns-1)*(num_sns)/2)])]);

%% covergence parameters
max_iter = 2000;
err_bnd = 1e-4;
err_max_ratio = 1-1e-10;
flag_converged = 0;
num_iter = 0;
use_convergence_criterion = 1;%0;%

%% Run algorithm
si_old = [0 0];
while err > err_bnd && flag_converged==0 && num_iter < max_iter
    
    if ~version_new
        [az_mat, el_mat, datas] = iter_ver0(az_mat, el_mat, datas, sns_pos, si_old, is2D, flag_sub_triang, flag_sub_plot);
    else
        [datas, SphMat, track_triang] = iter_ver1(datas, SphMat, sns_pos, si_old, is2D, flag_sub_triang, flag_sub_plot, track_triang);
    end
    
    %% errors - check for convergence
    errs = [];
    for si1 = 1:num_sns
        for si2 = si1+1:num_sns
            dif12_ORG2 = datas(si1).xyz' + repmat(sns_pos(:,si1)-sns_pos(:,si2),1,num_points) - datas(si2).xyz';
            erri = mean(sqrt(sum(dif12_ORG2.^2)));
            errs = [errs, erri];
        end
    end
    err_new = sum(errs)/((num_sns-1)*(num_sns)/2);
    if err_new/err > err_max_ratio
        flag_converged = use_convergence_criterion * 1; %stop iterations
    end
    err = err_new;
    %     disp(num2str([numel(err_vec), errs, err]));
    % do xyz
    % solve_wahba_svd
    % repeat
    err_vec = [err_vec err];
    
    %% errors vs. GT
    N_err_GT = numel(err_GT)+1;
    
    err_GT_i       = zeros(1,num_sns);
    err_GT_psi_i   = zeros(1,num_sns);
    err_GT_theta_i = zeros(1,num_sns);
    err_GT_phi_i   = zeros(1,num_sns);
    
    for si1 = 1:num_sns
        dif1_ORG0 = datas(si1).xyz' + repmat(sns_pos(:,si1),1,num_points) - GT_track;
        err_GT_i(si1) = mean(sqrt(sum(dif1_ORG0.^2)));
    end
    err_GT(N_err_GT).Cart = mean(err_GT_i);
    
    for si = 1:num_sns
        [~,C_new2old] = Wahba_XYZ(datas_original(si).xyz, datas(si).xyz, sns_pos(:,si), sns_pos(:,si));
        angles_sensor = dcm2angle(C_new2old)*pi/180;
        err_GT_psi_i(si)   = abs(angles_sensor(1)-GT_rot_psi(si)  );
        err_GT_theta_i(si) = abs(angles_sensor(2)-GT_rot_theta(si));
        err_GT_phi_i(si)   = abs(angles_sensor(3)-GT_rot_phi(si)  );
        
    end
    err_GT(N_err_GT).psi   = 1000*mean(err_GT_psi_i);
    err_GT(N_err_GT).theta = 1000*mean(err_GT_theta_i);
    err_GT(N_err_GT).phi   = 1000*mean(err_GT_phi_i);
    
    1;
    num_iter = num_iter + 1;
    
    if intermediate_plots && mod(num_iter,1)==0
        if plot_fig_100
            figure(100); hold on;
        end
        
        if exist('h','var')
            delete(h);
        end
        if version_new
            h = plott3(track_triang', '-','linewidth',2);
            disp(['iter: #',num2str(num_iter), ',  ', num2str(mean(sqrt(sum((GT_track - track_triang').^2)))) ]);
        else
            xyz_data = datas(1).xyz' + repmat(sns_pos(:,1),1,num_points);
            h = plott3(xyz_data', '-','linewidth',2);
            disp(['iter: #',num2str(num_iter), ',  ', num2str(mean(sqrt(sum((GT_track - xyz_data).^2)))) ]);
        end
        
        if isempty(M)
            M = getframe(gcf);
        else
            M(numel(M)+1) = getframe(gcf);
        end
        drawnow;
    end
    
    1;
end
writerObj = VideoWriter('E:\avi\myVideo.avi');
writerObj.FrameRate = 20;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
% for i=1:10
%     % convert the image to a frame
%     frame = M(1) ;
%     writeVideo(writerObj, frame);
% end
for i=1:length(M)
    % convert the image to a frame
    frame = M(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);



if plot_fig_100
    figure(100); hold on;
end
for si = 1:num_sns
    datas(si).xyz_bat(:,1) = datas(si).xyz(:,1) + sns_pos(1,si);
    datas(si).xyz_bat(:,2) = datas(si).xyz(:,2) + sns_pos(2,si);
    datas(si).xyz_bat(:,3) = datas(si).xyz(:,3) + sns_pos(3,si);
    
    
    if plot_fig_100
        plott3(datas(si).xyz_bat, '^-');
        drawnow;
    end
end
%--------------------------------------
% Compute the rotation matrices
angles = [];
for si = 1:num_sns
    [~,C_new2old] = Wahba_XYZ(datas_original(si).xyz, datas(si).xyz, sns_pos(:,si), sns_pos(:,si));
    angles_sensor = dcm2angle(C_new2old);
    angles = [angles, angles_sensor'];
end

end

%--------------------------------------------------------------------------
function disp_errs(datas,sns_pos,si1_in,si2_in,flag_plot)
num_sns   = length(sns_pos(1,:));
num_points = size(datas(1).xyz,1);
if flag_plot == 1
    errs = [];
    for si1 = 1:num_sns
        for si2 = si1+1:num_sns
            dif12_ORG2 = datas(si1).xyz' + repmat(sns_pos(:,si1)-sns_pos(:,si2),1,num_points) - datas(si2).xyz';
            erri = mean(sqrt(sum(dif12_ORG2.^2)));
            errs = [errs, erri];
        end
    end
    err = sum(errs)/((num_sns-1)*(num_sns)/2);
    disp(['     ', num2str(si1_in), '-->', num2str(si2_in), ':  ', num2str([errs, err])]);
end
end

%--------------------------------------------------------------------------
function [az_mat, el_mat, datas] = Data2Data_Triangulate_3sens(az_mat, el_mat, datas, sns_pos, fixed_sensor)
num_sns = length(sns_pos(1,:));

for si = 1:num_sns
    if ismember(si,fixed_sensor)
        [az_new,el_new,~] = cart2sph(datas(si).xyz(:,1),datas(si).xyz(:,2),datas(si).xyz(:,3));
        az_mat(:,si) = az_new;
        el_mat(:,si) = el_new;
    end
end

[datas, dops, SphMat, ~] = TriangulateManyGeneric(az_mat, el_mat, 0*az_mat, sns_pos, 1);
end

%--------------------------------------------------------------------------
function     [az_mat, el_mat, datas] = iter_ver0(az_mat, el_mat, datas, sns_pos, si_old, is2D, flag_sub_triang, flag_sub_plot)
num_sns   = length(sns_pos(1,:));
%%
for si1 = 1:num_sns
    for si2 = si1+1:num_sns
        
        if ismember(si1, si_old) %si1 is more updated. Now, first correct si2 to si1
            directions = [2 1 2 1 2 1];
        else
            directions = [1 2 1 2 1 2];
        end
        si_old = [si1 si2];
        %                 directions = [1 2 1 2 1 2];
        
        for iter_i = 1:2
            direction_i = directions(iter_i);
            switch direction_i
                case 1
                    %% si1 --> si2 
                    datas(si1).xyz = Wahba_XYZ(datas(si2).xyz, datas(si1).xyz, sns_pos(:,si2), sns_pos(:,si1)); %wahba
                    if is2D && flag_sub_triang==2
                        [az_mat, el_mat, datas] = Data2Data_Triangulate_3sens(az_mat, el_mat, datas, sns_pos, si1);
                    end
                    disp_errs(datas,sns_pos,si1,si2,flag_sub_plot);
                    
                case 2
                    %% si2 --> si1 
                    datas(si2).xyz = Wahba_XYZ(datas(si1).xyz, datas(si2).xyz, sns_pos(:,si1), sns_pos(:,si2));
                    if is2D && flag_sub_triang==1
                        [az_mat, el_mat, datas] = Data2Data_Triangulate_3sens(az_mat, el_mat, datas, sns_pos, [si1, si2]);
                    elseif is2D && flag_sub_triang==2
                        [az_mat, el_mat, datas] = Data2Data_Triangulate_3sens(az_mat, el_mat, datas, sns_pos, si2);
                    end
                    disp_errs(datas,sns_pos,si2,si1,flag_sub_plot);
            end
        end
        
        
    end
end

if is2D && flag_sub_triang==0
    [az_mat, el_mat, datas] = Data2Data_Triangulate_3sens(az_mat, el_mat, datas, sns_pos, [1:num_sns]);
end
end

%--------------------------------------------------------------------------
function     [datas, SphMat, track_triang] = iter_ver1(datas, SphMat, sns_pos, si_old, is2D, flag_sub_triang, flag_sub_plot, track_triang)
num_points = length(datas(1).xyz(:,1));
num_sns   = length(sns_pos(1,:));

%% stage 1 - 2Dcase - a. triangulate (az,el) samples from all sensors in order to get ragnes
%                     b. calc (x,y,z) from (az,el,range)
if is2D
    sensor_type = 1;
else
    sensor_type = 3;
end
[datas, dops, SphMat, track_triang] = TriangulateManyGeneric(SphMat.az_mat, SphMat.el_mat, SphMat.rng_mat, sns_pos, sensor_type);

%% stage 2 - calc rotation matrices using GT_hat
for si1 = 1:num_sns
    % GT_hat --> si1 
    [datas(si1).xyz] = Wahba_XYZ(track_triang, datas(si1).xyz, zeros(3,1), sns_pos(:,si1)); %wahba
end

%% stage 3 - recalc az&el (cart2sph)
% if is2D
for si=1:num_sns
    [SphMat.az_mat(:,si),SphMat.el_mat(:,si),SphMat.rng_mat(:,si)] = cart2sph(datas(si).xyz(:,1), datas(si).xyz(:,2), datas(si).xyz(:,3));
end
% end

%% stage 4 - recalc GT_hat using rotation matrices (track_triang)
si=1;
track_triang = datas(si).xyz + repmat(sns_pos(:,si),1,num_points)';
for si=2:num_sns
    track_triang = track_triang + datas(si).xyz + repmat(sns_pos(:,si),1,num_points)';
end
track_triang = track_triang/num_sns;

end

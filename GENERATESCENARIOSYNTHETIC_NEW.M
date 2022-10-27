function track = GenerateScenarioSynthetic_new(seg_num, init_states, omega_sign, ang_freq, omega_deg)

track = GenerateTraj(init_states,seg_num,omega_sign, ang_freq, omega_deg);
track = track(1:3,:);
end



%--------------------------------------------------------------------------

function states = GenerateTraj(init_state,seg_num,omega_sign, ang_freq, omega_deg)

segment_dur_str = 200;180; % sec
segment_dur_circ = ang_freq;180;225; % sec
% omega_deg = 1; % degree/sec
states = [];

for i=1:seg_num
    
    if mod(i,2)~=0
        states_segment = GenerateCV_model(init_state, segment_dur_str);
    else
        if mod(i,4)==0
            omega_deg = omega_sign * omega_deg;
        end
        states_segment = GenerateCT_model(init_state, segment_dur_circ, omega_deg);
    end
    states = [states, states_segment];
    init_state = states(:,end);
   
end
end
%--------------------------------------------------------------------------

function states_segment = GenerateCV_model(init_state, segment_dur)
states_segment = [init_state];
curr_state = init_state;
dt = 1;
A = [1 0 0 dt 0 0;...
     0 1 0 0 dt 0;...
     0 0 1 0 0 dt;...
     0 0 0 1 0 0;...
     0 0 0 0 1 0;...
     0 0 0 0 0 1];

for i=1:segment_dur
    
    new_state = A * curr_state;
    curr_state = new_state;
    states_segment = [states_segment, new_state];
    
end

end
%--------------------------------------------------------------------------

function states_segment = GenerateCT_model(init_state, segment_dur, omega_deg)
states_segment = [init_state];
curr_state = init_state;
dt = 1;
omega = omega_deg * pi/180; 
sin_omega = sin(omega*dt);
cos_omega = cos(omega*dt);

A = [1 0 0 sin_omega/omega -(1-cos_omega)/omega 0 ;...
     0 1 0 (1-cos_omega)/omega sin_omega/omega 0 ;...
     0 0 1 0 0 dt ;...
     0 0 0 cos_omega -sin_omega 0 ;...
     0 0 0 sin_omega cos_omega 0 ;...
     0 0 0 0 0 1];

for i=1:segment_dur
    
    new_state = A * curr_state;
    curr_state = new_state;
    states_segment = [states_segment, new_state];
    
end
end
%--------------------------------------------------------------------------

function track_state = GenerateTracksFromTraj(gt_states,err_bat,bias_bat,rot_bat)

num_sensors = length(err_bat);
track_gt = gt_states;
track_state = nan(size(track_gt));
trk_len = length(track_gt(1,:));
alpha = 0.85;
sigma = sqrt(err_bat^2/(1-alpha^2));

wp_now = zeros(3,1);
wv_now = zeros(3,1);

A = angle2dcm(rot_bat,rot_bat,rot_bat);

for i=1:trk_len
    
    wp_new = alpha * wp_now + sigma * randn(3,1) + bias_bat;
    wv_new = alpha * wv_now + sigma * randn(3,1);
    track_state(:,i) = [A*track_gt(1:3,i);track_gt(4:6,i)] + [wp_new;wv_new];
    wp_now = wp_new;
    wv_now = wv_new;
end

        
% figure
% plot(gt_states(4:6,:)')
% hold on; grid
% plot(track_state(4:6,:)')


end
%--------------------------------------------------------------------------

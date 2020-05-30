clc;
close all;
clear all;

%clearvars -except positions step_events beac_rssi_fixed_filtered beac_rssi_activity_filtered beac_motion trajectory_parts randn_x randn_y
participant = 1; %to choose the participant, from participant 1 to participant 8

[beac_rssi_fixed_filtered, beac_rssi_activity_filtered, beac_motion] = Load_data(participant);

load(step_events); % provide the line numbers in which Toe-off occur

step_events = [1, step_events, length(beac_motion)]; % Add '1' and the lenght of the dataset as points of support
positions=ZUPT_KF(acc, gyr_unbiased, gravity, acc_mean, fs, Step_events, participant, idx_fig+4);


% Get strides where the beacons are moved. That indicates that is the
% begining of a new trajectory
[trajectory_parts] = GetTrajectoryParts(beac_motion,step_events);


figure(participant)
xlim([2,18.8])
ylim([1.5, 9])
xlabel('x(meters)')
ylabel('y(meters)')
hold on
map
%grid on
distance_errors = zeros(10,10);

for i_iteration=1:10
disp(i_iteration)
N_PARTICLES = 600;
N_LM = 10;  % Number of landmarks(beacons)
LM_SIZE = 2; %landmark size (beacon's position in  two dimensions)
THRESHOLD_RESAMPLE = N_PARTICLES / 2;
STD_PERSON_POSITION = 0.1;

rssi_room = beac_rssi_fixed_filtered(:,1);
rssi_kitchen = beac_rssi_fixed_filtered(:,2);
rssi_bathroom = beac_rssi_fixed_filtered(:,3);
rssi_dining = beac_rssi_fixed_filtered(:,4);
rssi_living = beac_rssi_fixed_filtered(:,5);

rssi_door = beac_rssi_activity_filtered(:,1);
rssi_toilet = beac_rssi_activity_filtered(:,2);
rssi_broom = beac_rssi_activity_filtered(:,3);
rssi_pitcher = beac_rssi_activity_filtered(:,4);
rssi_brush = beac_rssi_activity_filtered(:,5);


%% MAIN LOOP
errors=[];
for i_trajectory_part=1:size(trajectory_parts,1)

    initial_step = trajectory_parts(i_trajectory_part,1);
    final_step = trajectory_parts(i_trajectory_part,2);
    pos_trajectory_part = positions(initial_step: final_step,:);
    
    if i_trajectory_part ==2
        x_initial = pos_trajectory_part(1, 1);
        y_initial = pos_trajectory_part(1, 2);
    else 
        x_initial = trajectory_parts(i_trajectory_part,3);
        y_initial = trajectory_parts(i_trajectory_part,4);
    end
    
    x_final =  trajectory_parts(i_trajectory_part,5);
    y_final =  trajectory_parts(i_trajectory_part,6);
        
    correction_x = pos_trajectory_part(1, 1) - x_initial;
    correction_y = pos_trajectory_part(1, 2) - y_initial;
    
    pos_trajectory_part(:, 1:2) = pos_trajectory_part(: ,1:2) - [correction_x correction_y];

    
    % CREATION OF PARTICLES
    particles = [];
    trajectories = cell(1,N_PARTICLES);
    landmarks_figures = cell(1,N_LM);
    particles_weights = ones(1,N_PARTICLES) / N_PARTICLES;
    for i_particle=1:N_PARTICLES
        particles = [particles, Particle(N_PARTICLES,N_LM,LM_SIZE, x_initial, y_initial)];
        trajectory = plot(nan, nan, 'color', [.9 .9 .9]);
        trajectories{i_particle} = trajectory;
    end
    for i_landmark=1:N_LM
        landmark_figure = plot(nan, nan,'d');
        landmarks_figures{i_landmark} = landmark_figure;
    end
    
    particles = init_beacons_position(particles,N_PARTICLES);
    current_best_particle = N_PARTICLES;
    previous_best_particle = N_PARTICLES;
    
    t_x=[];
    t_y=[];
    
    for i_step=2:size(pos_trajectory_part,1)

        % PREDICTION
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = pos_trajectory_part(i_step, 1);
        y = pos_trajectory_part(i_step, 2);
        for i_particle=1:N_PARTICLES
            
            dist_x = (x - x_initial) + randn() * STD_PERSON_POSITION;
            dist_y = (y - y_initial) + randn() * STD_PERSON_POSITION;
            particles(i_particle).X = particles(i_particle).X + dist_x;
            particles(i_particle).Y = particles(i_particle).Y + dist_y;
            
            particles(i_particle).T_x = [particles(1,i_particle).T_x, particles(i_particle).X];
            particles(i_particle).T_y = [particles(1,i_particle).T_y, particles(i_particle).Y];
            
            trajectory_x_data = [particles(i_particle).T_x];
            trajectory_y_data = [particles(i_particle).T_y];
            
            if i_trajectory_part==1
                t_x = [t_x x] ;
                t_y = [t_y y] ;
                set(trajectories{i_particle},'XData',t_x, 'YData',t_y, 'color', 'b');
            else
                set(trajectories{i_particle},'XData',trajectory_x_data, 'YData',trajectory_y_data, 'color', [.9 .9 .9]);
            end
            
        end


        
        x_initial = x;
        y_initial = y;
        pause(0.01)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % UPDATE
        initial_sample = step_events(i_step + initial_step -1);
        final_sample = step_events(i_step + initial_step);
        % Do not include RSSI measurements while person is still
        if (final_sample - initial_sample) > 1024
            initial_sample = final_sample - 1024;
        end
        
        rssis = [rssi_room(initial_sample:final_sample), rssi_kitchen(initial_sample:final_sample), rssi_bathroom(initial_sample:final_sample),...
                 rssi_dining(initial_sample:final_sample), rssi_living(initial_sample:final_sample)...
                 rssi_door(initial_sample:final_sample), rssi_toilet(initial_sample:final_sample), rssi_broom(initial_sample:final_sample),...
                 rssi_pitcher(initial_sample:final_sample), rssi_brush(initial_sample:final_sample)];
        
        for i_rssis=1:size(rssis,1)
            non_zero_indexes = find(rssis(i_rssis,:));
            if ~isempty(non_zero_indexes)
                
                for j_rssi=non_zero_indexes
                    rssi = rssis(i_rssis,j_rssi);
                    current_sample = (initial_sample -1 +i_rssis);
                    
                    if j_rssi <= 5 %Stationary beacons
                        if rssi > -85
                            [particles, particles_weights, trajectories, current_best_particle] = ...
                                update_resample(particles, rssi, j_rssi, particles_weights, N_PARTICLES, trajectories, THRESHOLD_RESAMPLE, previous_best_particle);
                        end
                    end                   
                end
            end
        end        
    end
    
    

    for i_particle=1:N_PARTICLES
        if i_particle ~= current_best_particle
        set(trajectories{i_particle},'XData',[], 'YData',[]);
        end
    end
    set(trajectories{current_best_particle},'color', 'b');
    
    error = sqrt((particles(current_best_particle).X - x_final)^2 + (particles(current_best_particle).Y - y_final)^2);
    errors = [errors error];
    
%     % Plot the landmarks of the best particle.
%     for i_lm=1:N_LM
%         hold on
%         set(landmarks_figures{i_lm},'XData',particles(current_best_particle).Lm(i_lm,1), 'YData',particles(current_best_particle).Lm(i_lm,2));
%     end
    
end
distance_errors(i_iteration,1:length(errors)) = errors;
end




function [particles, particles_weights, trajectories, current_best_particle] = ...
    update_resample(particles, rssi, j_rssi, particles_weights, N_PARTICLES, trajectories, THRESHOLD_RESAMPLE, previous_best_particle)
    
    % To update
    [particles, particles_weights, current_best_particle] = update(particles, rssi, j_rssi, particles_weights, N_PARTICLES, previous_best_particle);

    % To resample
    neff = 1 / sum(sqrt(particles_weights));
    if neff < THRESHOLD_RESAMPLE
        %indexes = double(py.mifunc.stratified_resample(particles_weights))+1;
        indexes = double(py.mifunc.systematic_resample(particles_weights))+1;
        [particles, particles_weights, trajectories] = resample_from_index(particles, particles_weights, indexes, N_PARTICLES, trajectories);
        %particles_weights = ones(1,N_PARTICLES) / N_PARTICLES; % PARTICLES WEIGHTS ARE INITILIZED TO AND EQUAL VALUE

    end
end

function [particles, particles_weights, current_best_particle] = update(particles, rssi, lm_id, particles_weights, N_PARTICLES, previous_best_particle)
    
    % different models of beacons have different propagation parameters
    if lm_id == 5 % living
        distance_rssi = 10 ^ ((rssi + 78) / -10 * 0.6);  % Calculate distance to beacon
        R = distance_rssi / 10;
    else
        distance_rssi = 10 ^ ((rssi + 83) / -10 * 0.6);  % Calculate distance to beacon
        R = distance_rssi / 10;
    end

    for i_particle=1:N_PARTICLES
        % Get the distance between the person and the landmark
        dx = particles(i_particle).X - particles(i_particle).Lm(lm_id,1);
        dy = particles(i_particle).Y - particles(i_particle).Lm(lm_id,2);
        dist = sqrt((dx * dx) + (dy * dy));
        residual = distance_rssi - dist;

        % Compute Jacobians
        H = [-dx / dist, -dy / dist];

        % Compute covariance of the residual
        % covV = H * Cov_s * H^T + error
        HxCov = [particles(i_particle).LmP{1,lm_id}(1, 1) * H(1) + particles(i_particle).LmP{1,lm_id}(1, 2) * H(2),...
                 particles(i_particle).LmP{1,lm_id}(2, 1) * H(1) + particles(i_particle).LmP{1,lm_id}(2, 2) * H(2)];

        covV = (HxCov(1) * H(1)) + (HxCov(2) * H(2)) + R;

        % Calculate Kalman gain
        K_gain = [HxCov(1) * (1 / covV), HxCov(2) * (1.0 / covV)];

        % Calculate the new landmark position
        lm_x = particles(i_particle).Lm(lm_id,1) + (K_gain(1) * residual);
        lm_y = particles(1,i_particle).Lm(lm_id,2) + (K_gain(2) * residual);

        % Calculate the new covariance matrix of the landmark
        % cov_t = cov_t-1 - K * covV * K^T
        lm_P_aux = [[K_gain(1) * K_gain(1) * covV, K_gain(1) * K_gain(2) * covV];...
                    [K_gain(2) * K_gain(1) * covV, K_gain(2) * K_gain(2) * covV]];

        lm_P = [[particles(i_particle).LmP{1,lm_id}(1, 1) - lm_P_aux(1,1),...
                 particles(i_particle).LmP{1,lm_id}(1, 2) - lm_P_aux(1,2)];...
                [particles(i_particle).LmP{1,lm_id}(2, 1) - lm_P_aux(2,1),...
                 particles(i_particle).LmP{1,lm_id}(2, 2) - lm_P_aux(2,2)]];

        % Update landmark in particle
        particles(i_particle).Lm(lm_id, 1) = lm_x;
        particles(i_particle).Lm(lm_id, 2) = lm_y;
        particles(i_particle).LmP{1,lm_id} = lm_P;

        % Update particles weight
        particles_weights(i_particle) = py.mifunc.calculate_weight(particles_weights(i_particle), dist, covV, distance_rssi);

        %particles_weights(i_particle) = (particles_weights(i_particle) * ...
        %                                (1 / (covV * sqrt(2*pi))) * exp((-(distance_rssi - dist)^2) / (2 * covV^2))) + 1*exp(-300);
    end
    
    particles_weights = particles_weights / sum(particles_weights);
    best_particles = [];
    biggest_weight = -1;
    for i_particle=1:N_PARTICLES
        particles(i_particle).W = particles_weights(i_particle);
        
        % Get the best particle (particle with the biggest weight)
        if particles_weights(i_particle) >= biggest_weight
            best_particles = [best_particles, i_particle];
        end
        if find(best_particles==previous_best_particle)
            current_best_particle = previous_best_particle;
        else
            current_best_particle = best_particles(1); % for reproducibility the first item is always chosen.
        end
        
    end
end

function [new_particles, new_particles_weights, trajectories] = resample_from_index(particles, particles_weights, indexes, N_PARTICLES, trajectories)
    new_particles_weights = zeros(1,N_PARTICLES);
    new_particles = particles;
    for i_particle = 1:N_PARTICLES
        new_particles_weights(i_particle) = particles_weights(indexes(i_particle));
        new_particles(i_particle) = particles(indexes(i_particle));
    end
    
    new_particles_weights = new_particles_weights / sum(new_particles_weights);
    for i_particle=1:N_PARTICLES
        %new_particles(1,i_particle).W =
        %new_particles_weights(i_particle);%%        
        new_particles(i_particle).W = 1 / N_PARTICLES;
        new_particles_weights(i_particle) = 1 / N_PARTICLES;%%
        
        trajectory_x_data = [new_particles(i_particle).T_x];
        trajectory_y_data = [new_particles(i_particle).T_y];
        set(trajectories{i_particle},'XData',trajectory_x_data, 'YData',trajectory_y_data);

    end
    %pause(0.005)
end
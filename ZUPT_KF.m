function [positions]=ZUPT_KF(acc, gyr_unbiased, g, acc_mean, fs, Step_events, participant, idx_fig)
%% Read data from file.
% Data should include timestamps (seconds), 3 axis accelerations (m/s^2), 3
% axis gyroscopic rates of turn (rad/s).


data_size = size(acc,1);
Acc = acc';
Gyro = gyr_unbiased';


%% Initialise parameters.
% Orientation from accelerometers. Sensor is assumed to be stationary.
pitch =-atan2(acc_mean(1),sqrt(acc_mean(2)^2+acc_mean(3)^2));
roll = atan2(acc_mean(2),acc_mean(3));  
yaw = 0;

C = [cos(pitch)*cos(yaw) (sin(roll)*sin(pitch)*cos(yaw))-(cos(roll)*sin(yaw)) (cos(roll)*sin(pitch)*cos(yaw))+(sin(roll)*sin(yaw));
    cos(pitch)*sin(yaw)  (sin(roll)*sin(pitch)*sin(yaw))+(cos(roll)*cos(yaw))  (cos(roll)*sin(pitch)*sin(yaw))-(sin(roll)*cos(yaw));
    -sin(pitch) sin(roll)*cos(pitch) cos(roll)*cos(pitch)];
C_prev = C;

% Preallocate storage for heading estimate. Different from direction of
% travel, the heading indicates the direction that the sensor, and therefore
% the pedestrian, is facing.
heading = nan(1, data_size);
heading(1) = yaw;

% Gyroscope bias, to be determined for each sensor.
%  -- Defined above so we don't forget to change for each dataset. --

% Preallocate storage for accelerations in navigation frame.
acc_n = nan(3, data_size);
acc_n(:,1) = C*Acc(:,1);


% Preallocate storage for velocity (in navigation frame).
% Initial velocity assumed to be zero.
vel_n = nan(3, data_size);
vel_n(:,1) = [0 0 0]';

% Preallocate storage for position (in navigation frame).
% Initial position arbitrarily set to the origin.
pos_n = nan(3, data_size);
pos_n(:,1) = [0   0   0]';

% Preallocate storage for distance travelled used for altitude plots.
distance = nan(1,data_size-1);
distance(1) = 0;


% Error covariance matrix.
P = zeros(9);

% Process noise parameter, gyroscope and accelerometer noise.
sigma_omega = 0.02;
sigma_a = 0.03;

% ZUPT measurement matrix.
H = [zeros(3) zeros(3) eye(3)];

% ZUPT measurement noise covariance matrix.
sigma_v = 0.005;
R = diag([sigma_v sigma_v sigma_v]).^2;

% Gyroscope stance phase detection threshold.
gyro_threshold = 0.6;

%% Main Loop
for t = 2:data_size
    %%% Start INS (transformation, double integration) %%%
    dt = 1/fs;
    
    % Remove bias from gyro measurements.
    gyro_s1 = Gyro(:,t); %- gyro_bias;
    
    % Skew-symmetric matrix for angular rates
    ang_rate_matrix = [0   -gyro_s1(3)   gyro_s1(2);
        gyro_s1(3)  0   -gyro_s1(1);
        -gyro_s1(2)  gyro_s1(1)  0];
    
    % orientation estimation
    C = C_prev*(2*eye(3)+(ang_rate_matrix*dt))/(2*eye(3)-(ang_rate_matrix*dt));
    
    % Transforming the acceleration from sensor frame to navigation frame.
    acc_n(:,t) = 0.5*(C + C_prev)*Acc(:,t);
    
    % Velocity and position estimation using trapeze integration.
    vel_n(:,t) = vel_n(:,t-1) + ((acc_n(:,t) - [0; 0; g] )+(acc_n(:,t-1) - [0; 0; g]))*dt/2;
    pos_n(:,t) = pos_n(:,t-1) + (vel_n(:,t) + vel_n(:,t-1))*dt/2;
    
    % Skew-symmetric cross-product operator matrix formed from the n-frame accelerations.
    S = [0  -acc_n(3,t)  acc_n(2,t);
        acc_n(3,t)  0  -acc_n(1,t);
        -acc_n(2,t) acc_n(1,t) 0];
    
    % State transition matrix.
    F = [eye(3)  zeros(3,3)    zeros(3,3);
        zeros(3,3)   eye(3)  dt*eye(3);
        -dt*S  zeros(3,3)    eye(3) ];
    
    % Compute the process noise covariance Q.
    Q = diag([sigma_omega sigma_omega sigma_omega 0 0 0 sigma_a sigma_a sigma_a]*dt).^2;
    
    % Propagate the error covariance matrix.
    P = F*P*F' + Q;
    %%% End INS %%%
    
    % Stance phase detection and zero-velocity updates.
    if norm(Gyro(:,t)) < gyro_threshold
        %%% Start Kalman filter zero-velocity update %%%
        % Kalman gain.
        K = (P*(H)')/((H)*P*(H)' + R);
        
        % Update the filter state.
        delta_x = K*vel_n(:,t);
        
        % Update the error covariance matrix.
        P = (eye(9) - K*H)*P;
        
        % Extract errors from the KF state.
        attitude_error = delta_x(1:3);
        pos_error = delta_x(4:6);
        vel_error = delta_x(7:9);
        %%% End Kalman filter zero-velocity update %%%
        
        %%% Apply corrections to INS estimates. %%%
        % Skew-symmetric matrix for small angles to correct orientation.
        ang_matrix = -[0   -attitude_error(3,1)   attitude_error(2,1);
            attitude_error(3,1)  0   -attitude_error(1,1);
            -attitude_error(2,1)  attitude_error(1,1)  0];
        
        % Correct orientation.
        C = (2*eye(3)+(ang_matrix))/(2*eye(3)-(ang_matrix))*C;
        
        % Correct position and velocity based on Kalman error estimates.
        vel_n(:,t)=vel_n(:,t)-vel_error;
        pos_n(:,t)=pos_n(:,t)-pos_error;
    end
    heading(t) = atan2(C(2,1), C(1,1)); % Estimate and save the yaw of the sensor (different from the direction of travel). Unused here but potentially useful for orienting a GUI correctly.
    C_prev = C; % Save orientation estimate, required at start of main loop.
    
    % Compute horizontal distance.
    distance(1,t) = distance(1,t-1) + sqrt((pos_n(1,t)-pos_n(1,t-1))^2 + (pos_n(2,t)-pos_n(2,t-1))^2);
end



%Rotate position estimates and plot.
angles = [-65 -130 -75 -75 -75 -65 0 -75 -50 -165 -80]; % Rotation angle required to achieve an aesthetic alignment of the figure.
angle = angles(participant);
rotation_matrix = [cosd(angle) -sind(angle);
    sind(angle) cosd(angle)];
pos_r = zeros(2,length(pos_n));
for idx = 1:length(pos_n)
    pos_r(:,idx) = rotation_matrix*[pos_n(1,idx) pos_n(2,idx)]';
end
pos_r(1,:) = pos_r(1,:) + 3;
pos_r(2,:) = pos_r(2,:) + 7.5;
pos_r(3,:) = pos_n(3,:)';
positions_kf = [pos_r(1,:); pos_r(2,:); pos_r(3,:)]'; % minus to pass from west to east

% Plot in step basis
number_of_steps = length(Step_events);
positions = zeros(number_of_steps+1,3);
positions(1,:) = [3 7.5 0]; % initial position
distances = zeros(number_of_steps,1);
for i_step=1:number_of_steps
    sample_step = Step_events(i_step);
    positions(i_step+1,:) = [positions_kf(sample_step,1), positions_kf(sample_step,2), positions_kf(sample_step,3)]; % minus to pass from west to east
    distances(i_step,1) = distance(sample_step);
end
end
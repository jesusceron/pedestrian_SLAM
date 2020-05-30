function [beac_rssi_fixed_filtered, beac_rssi_activity_filtered, beac_motion] = Load_data(participant)


    %% Load files    
    
    Acc_Gyr_Beac_data= readtable(strcat('\dataset_synchronized\',int2str(participant),'.csv'));

    gyro_bias_i_c_s = [820 5500 4000 1 5000 3600 16120 1000 500 4000 1];
    gyro_bias_f_c_s = [1420 7500 8000 7500 6300 4300 17220 3000 2000 5000 3000];
    start_data_sample = [8080 8700 9500 8600 7000 5300 18800 9150 5600 5880 3900];
    final_data_sample = [63930 60200 88900 71200 72400 55450 114800 72600 70000 54630 62720];

    gyro_i_c_s = gyro_bias_i_c_s(participant); % initial calibration sample
    gyro_f_c_s = gyro_bias_f_c_s(participant); % initial calibration sample

    start_sample = start_data_sample(participant);
    final_sample = final_data_sample(participant);
    fs = 204.8; % IMU sample rate

    %% Load data.

    % Load IMU data
    acc_complete = table2array(Acc_Gyr_Beac_data(:,23:25));
    gyro_complete = table2array(Acc_Gyr_Beac_data(:,26:28));

    % Load BEACONS data.
    beac_rssi_fixed = table2array(Acc_Gyr_Beac_data(:,...
        {'RSSIs_beacon_1',...   % Coconut:  Room
        'RSSIs_beacon_2',...    % Mint:     Kitchen
        'RSSIs_beacon_3',...    % Ice:      Bathroom
        'RSSIs_beacon_4',...    % Blueber:  Dining room
        'RSSIs_beacon_5'}));    % P2:       Living room

    beac_rssi_activity = table2array(Acc_Gyr_Beac_data(:,...
        {'RSSIs_beacon_6',...   % P1:       Door
        'RSSIs_beacon_7',...    % B2:       Toilet lid
        'RSSIs_beacon_8',...    % B1:       Broom
        'RSSIs_beacon_9',...    % G2:       Pitcher
        'RSSIs_beacon_10'}));   % G1:       Hair brush

    beac_motion = table2array(Acc_Gyr_Beac_data(:,{'motion_state_beacon_6',...
        'motion_state_beacon_7',...
        'motion_state_beacon_8',...
        'motion_state_beacon_9',...
        'motion_state_beacon_10'}));

    %% Data pre-processing
    w = gyro_i_c_s : gyro_f_c_s;

    acc = zeros(length(acc_complete),3);
    acc(:,1) = acc_complete(:,2); 
    acc(:,2) = acc_complete(:,3); % y = z
    acc(:,3) = acc_complete(:,1); 

    gyr = zeros(length(gyro_complete),3);
    gyr(:,1) = gyro_complete(:,2); 
    gyr(:,2) = gyro_complete(:,3); 
    gyr(:,3) = gyro_complete(:,1); 
    gyr = deg2rad(gyr);

    gravity=mean(sqrt(acc(w,1).^2 + acc(w,2).^2 + acc(w,3).^2));  % m/s^2
    acc_mean=[mean(acc(w,1)),mean(acc(w,2)),mean(acc(w,3))];

    % Remove bias Gyro
    bias_gyr=[mean(gyr(w, 1)), mean(gyr(w, 2)), mean(gyr(w, 3))];
    gyr_unbiased=[gyr(:,1) - bias_gyr(1),...
        gyr(:,2) - bias_gyr(2),...
        gyr(:,3) - bias_gyr(3)];

    acc = acc(start_sample:final_sample,:);
    gyr_unbiased = gyr_unbiased(start_sample:final_sample,:);
    beac_rssi_fixed = beac_rssi_fixed(start_sample:final_sample,:);
    beac_rssi_activity = beac_rssi_activity(start_sample:final_sample,:);
    beac_motion = beac_motion(start_sample:final_sample,:);

    % Filtering the beacons data in Python
    [beac_rssi_fixed_filtered, beac_rssi_activity_filtered] = rssiKF(beac_rssi_fixed,beac_rssi_activity);

end

    
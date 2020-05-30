function [trajectory_parts] = GetTrajectoryParts(beac_motion,step_events)
% Get strides where the beacons are moved. That indicates that is the
% begining of a new trajectory

gap = 508; % movement readings have a delay of about 2 seconds (204.8Hz = 508 samples)

door_moving = [beac_motion(gap:end,1); zeros(gap-1,1)];
toilet_moving = [beac_motion(gap:end,2); zeros(gap-1,1)];
broom_moving = [beac_motion(gap:end,3); zeros(gap-1,1)];
pitcher_moving = [beac_motion(gap:end,4); zeros(gap-1,1)];
brush_moving = [beac_motion(gap:end,5); zeros(gap-1,1)];

for i_beacon=1:5
    
    switch i_beacon
        case 1
            door_movement = find(door_moving);
            sample_first_door_movement = door_movement(1);
            for i_door_movement=1:size(door_movement,1)
                if door_movement(i_door_movement) > sample_first_door_movement + 20000 % % add 2000; gap between the event of entering and leaving the apartment
                    sample_second_door_movement = door_movement(i_door_movement);
                end
            end
            
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_first_door_movement
                    first_step_in_door = i_step;
                    door_position = [5 5.5];
                    break
                end
            end
            
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_second_door_movement
                    second_step_in_door = i_step;
                    door_position = [6.2 6];
                    break
                end
            end
            
        case 2
            toilet_movement = find(toilet_moving);
            sample_toilet_movement = toilet_movement(1);
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_toilet_movement
                    step_in_toilet = i_step;
                    toilet_position = [8.8 8.2];
                    break
                end
            end
        case 3
            broom_movement = find(broom_moving);
            sample_broom_movement = broom_movement(1);
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_broom_movement
                    step_in_broom = i_step ;
                    broom_position = [7.5 8];
                    break
                end
            end
        case 4
            pitcher_movement = find(pitcher_moving);
            sample_pitcher_movement = pitcher_movement(1);
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_pitcher_movement
                    step_in_pitcher = i_step;
                    pitcher_position = [12.6 7.4];
                    break
                end
            end
        case 5
            brush_movement = find(brush_moving);
            sample_brush_movement = brush_movement(1);
            for i_step=1:length(step_events)
                if step_events(i_step) > sample_brush_movement
                    step_in_brush = i_step;
                    brush_position = [8.7 6.8];
                    break
                end
            end
    end
end

a = [first_step_in_door; second_step_in_door; step_in_toilet; step_in_broom; step_in_pitcher; step_in_brush];
pos = [door_position; door_position; toilet_position; broom_position; pitcher_position; brush_position];
[in, b] = sort(a);
p = [in pos(b,:)];

init_step = [1; p(:,1)];
final_step = [p(:,1);length(step_events)-1];
init_pos = [[3, 7.5]; p(:,2:3)];
final_pos = [p(:,2:3);[3, 7.5]];

trajectory_parts = [init_step, final_step, init_pos, final_pos];
function particles = init_beacons_position(particles,N_PARTICLES)
% Initialization of beacons position
%   at the very beginning the position of the beacons is unknown.
%   the covariance matrix is initialized with big variance values.

% Initialize Lm[1] (beacon in the room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(1,:) = [0, 0];
        particles(1,i_particle).LmP{1} = [[100, 0]; [0, 100]];
    end

    % Initialize Lm[2] (beacon in the kitchen)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(2,:) = [0, 0];
        particles(1,i_particle).LmP{2} = [[100, 0]; [0, 100]];
    end

    % Initialize Lm[3] (beacon in the bathroom)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(3,:) = [0, 0];
        particles(1,i_particle).LmP{3} = [[100, 0]; [0, 100]];
    end

    % Initialize Lm[4] (beacon in the dining room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(4,:) = [0, 0];
        particles(1,i_particle).LmP{4} = [[100, 0]; [0, 100]];
    end

    % Initialize Lm[5] (beacon in the living room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(5,:) = [0, 0];
        particles(1,i_particle).LmP{5} = [[100, 0]; [0, 100]];
    end

    % Initialize Lm[6] (beacon in the door)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(6,:) = [0, 0];
        particles(1,i_particle).LmP{6}= [[100, 0]; [0, 100]];
    end

    % Initialize Lm[7] (beacon in the toilet)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(7,:) = [0, 0];
        particles(1,i_particle).LmP{7} = [[100, 0]; [0, 100]];
    end

    % Initialize Lm[8] (beacon in the broom)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(8,:) = [0, 0];
        particles(1,i_particle).LmP{8} =  [[100, 0]; [0, 100]];
    end

    % Initialize Lm[9] (beacon in the pitcher)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(9,:) = [0, 0];
        particles(1,i_particle).LmP{9} =  [[100, 0]; [0, 100]];
    end

    % Initialize Lm[10] (beacon in the brush)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(10,:) = [0, 0];
        particles(1,i_particle).LmP{10} =  [[100, 0]; [0, 100]];
    end

end
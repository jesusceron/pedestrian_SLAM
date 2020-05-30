classdef Particle
    %Particle. Trajectory particles
    %   Each particle represents a different system state.
    %   
    
    properties
        W
        X
        Y
        Lm
        LmP
        T_x
        T_y
    end
    
    methods
        function obj = Particle(N_PARTICLES,N_LM,LM_SIZE,x,y)
            
            obj.W = 1/N_PARTICLES;
            obj.X = x; %position in x
            obj.Y = y; %position in y
            obj.Lm = zeros(N_LM, LM_SIZE); %set of the positions of the beacons for this particle
            obj.LmP = {}; % Covarian matrix
            for i=1:N_LM
                obj.LmP{i} = zeros(LM_SIZE, LM_SIZE);
            end
            % historical of positions (person's trajectory)
            obj.T_x = x; 
            obj.T_y = y;            
            
        end

    end
end


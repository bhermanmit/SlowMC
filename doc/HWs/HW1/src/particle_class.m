classdef particle_class
    %PARTICLE_CLASS defines a particle type for slowing down
    %   this class defines a particle along with its methods that will
    %   be used in the Monte Carlo slowing down code
    
    properties (SetAccess = private, GetAccess = public)
        E
        alive
    end
    
    methods
        
        % Constructor
        function obj=particle_class()    
            obj.alive = 1;            
        end
        
        % kill particle
        function obj=kill(obj)
            obj.alive = 0;
        end
        
        % sets particle energy
        function obj = set_energy(obj,E)
            obj.E = E;
        end
        
    end
    
end


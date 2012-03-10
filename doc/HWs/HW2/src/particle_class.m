classdef particle_class
    %PARTICLE_CLASS defines a particle type for slowing down
    %   this class defines a particle along with its methods that will
    %   be used in the Monte Carlo slowing down code
    
    properties (SetAccess = private, GetAccess = public)
        E
        alive
        n_coll_therm
        n_coll
        n_abs
        time
        mass            % neutron mass in kg
    end
    
    methods
        
        % Constructor
        function obj=particle_class()    
            obj.alive = 1;
            obj.n_coll_therm = 0;
            obj.n_coll = 0;
            obj.n_abs = 0;
            obj.time = 0;
            obj.mass = 1.674927351e-27;
        end
        
        % kill particle
        function obj=kill(obj)
            obj.alive = 0;
        end
        
        % sets particle energy
        function obj = set_energy(obj,E)
            obj.E = E;
        end
        
        % accumulate neutron lifetime
        function obj = add_time(obj,xs)
           
            % calculate velocity
            v = sqrt(2*obj.E*1.60217653e-13/obj.mass);
            
            % accumulate lifetime
            obj.time = obj.time + 1/(v*xs);
            
        end
        
        % increment the collision vars
        function obj = collision(obj,reaction)
           
            % bank collision info
            if obj.E > 1e-6
               obj.n_coll_therm = obj.n_coll_therm + 1; 
            end
            obj.n_coll = obj.n_coll + 1;
            
            % bank absorption if absorbed
            if strcmp(reaction,'capture')
                obj.n_abs = obj.n_abs + 1;
            end
        end
        
        % increment the abs var
        function obj = increment_abs(obj)
           obj.n_abs = obj.n_abs + 1; 
        end
        
    end
    
end
classdef cross_section_class
    %CROSS_SECTION keeps track of ENDF xs data sets
    
    properties
        egrid
        xs
        name
        N
        totxs
    end
    
    methods
        
        % constructor
        function obj = cross_section_class(name,N)
           
            % get name
            obj.name = name;
            
            % load xs data file expecting E and xs
            disp(strcat('Loading xsdata file: ',name,'.m'));
            load(horzcat('./xsdata/',name));
            
            % set energy grid and xs
            obj.egrid = E;
            obj.xs = xs;
            
            % set number density
            obj.N = N;
            
            % save macroscopic total cross section
            obj.totxs = N*sum(obj.xs,2);
            
        end
    end
    
end


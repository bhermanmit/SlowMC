classdef cross_section_class
    %CROSS_SECTION keeps track of ENDF xs data sets
    
    properties
        egrid
        xs
        name
    end
    
    methods
        
        % constructor
        function obj = cross_section_class(name)
           
            % get name
            obj.name = name;
            
            % load xs data file expecting E and xs
            disp(strcat('Loading xsdata file: ',name,'.m'));
            load(name);
            
            % set energy grid and xs
            obj.egrid = E;
            obj.xs = xs;
            
        end
    end
    
end


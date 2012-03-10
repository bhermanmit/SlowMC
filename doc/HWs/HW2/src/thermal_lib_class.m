classdef thermal_lib_class
    %THERMAL_LIB_CLASS defines the thermal scattering library object
    
    properties
        Erat
        cdf
        name
        kT
    end
    
    methods
        
        % constructor
        function obj = thermal_lib_class(name)
           
            % get name
            obj.name = name;
            
            % load xs data file expecting E and xs
            disp(horzcat('Loading thermal cdf file: ',name,'_therm'));
            load(horzcat('./pdfdata/',name,'_therm'));
            
            % set energy grid and xs
            obj.Erat = Erat;
            obj.cdf = thermalcdf;
            obj.kT = kT;
            
        end
    end
    
end
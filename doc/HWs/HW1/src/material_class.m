classdef material_class
    %MATERIAL_CLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        A
        alpha
        n_isotopes
        isotopes
        iso_list
    end
    
    methods
        
        % constructor
        function obj = material_class(A)
            
            % set vars
            obj.A = A;
            obj.alpha = ((A-1)/(A+1))^2;
            obj.n_isotopes = 0;
            
        end
        
        % import xsdata file
        function obj = load_isotope(obj,name)
            
            % increment number of isotopes
            obj.n_isotopes = obj.n_isotopes + 1;
            
            % call constructor of xsdata
            obj.isotopes{obj.n_isotopes} = cross_section_class(name);
            
            % append to list
            obj.iso_list{obj.n_isotopes} = name;
            
        end
        
    end
    
end


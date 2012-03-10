classdef material_class
    %MATERIAL_CLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        A
        alpha
        n_isotopes
        isotopes
        iso_list
        totxs
        xs_ref = 7;
        E_ref = 0.025e-6;
    end
    
    methods
        
        % constructor
        function obj = material_class(A)
            
            % set vars
            obj.A = A;
            obj.alpha = ((A-1)/(A+1))^2;
            obj.n_isotopes = 0;
            obj.totxs = 0;
            
        end
        
        % import xsdata file
        function obj = load_isotope(obj,name,N)
            
            % increment number of isotopes
            obj.n_isotopes = obj.n_isotopes + 1;
            
            % call constructor of xsdata
            obj.isotopes{obj.n_isotopes} = cross_section_class(name,N);
            
            % append to list
            obj.iso_list{obj.n_isotopes} = name;
            
        end
        
        % get total xs
        function obj = compute_macro_totxs(obj,E)
            
            % reset totxs
            obj.totxs = 0;
            
            % begin loop around isotopes
            for i = 1:obj.n_isotopes
                
                % get macro total cross section
                xs_tot = interp1(obj.isotopes{i}.egrid,obj.isotopes{i}.xs,E,'linear','extrap');
            
                % get 1/v absorption cross section
                xs_a = sqrt(obj.E_ref/E)*obj.xs_ref;
            
                % compute total xs
                obj.totxs = obj.totxs + xs_tot + xs_a;
                
            end
            
            
        end
        
    end
    
end


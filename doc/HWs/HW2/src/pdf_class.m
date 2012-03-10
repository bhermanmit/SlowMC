classdef pdf_class
    
    % PDF_CLASS contains information about the distribution function
    %   contains all of the random number sampling and pdfs needed by MC
    
    properties (SetAccess = private, GetAccess = public)
      
        rng    % random number generator object
        E      % new sampled energy
        egrid  % energy grid for chi pdf
        chicdf % cdf for chi
        n_thermal % number of thermal scattering libs
        thermal_list % list of thermal isotopes
        thermallibs % array of thermal scattering objs
        kT = 8.6173324e-5*300*10^-6; % boltzmann constant times temperature
        xs_ref = 7;
        E_ref = 0.025e-6;
        reaction
        
    end
    
    methods
        
        % constructor
        function obj = pdf_class(seed)
           
            % initialize random number generator with seed
            obj.rng = RandStream('mcg16807','Seed',seed);
            
            % create cdf
            obj = obj.watt_fission_cdf();
            
            % set thermal iso to zero
            obj.n_thermal = 0;
            
        end
        
        % sample energy
        function obj = sample_collision_energy(obj,E,alpha)
           
            if E > 4e-6
                obj.E = E - E*(1-alpha)*obj.rng.rand;
            else
                
                % get a random number
                rn = obj.rng.rand;
                
                % create vector of energies
                Evec = zeros(1,size(obj.thermallibs{1}.cdf,2));
                
                % loop around a thermal kernels and get indices
                for i = 1:length(Evec)
                    
                    % get index
                    idx = find(obj.thermallibs{1}.cdf(:,i).* ...
                         (obj.thermallibs{1}.cdf(:,i) > rn),1,'first');
                    
                    % get possible outgoing energy ratio
                    Evec(i) = obj.thermallibs{1}.Erat(idx);
                    
                end
                
                % interpolate energy ration to get correct outgoing energy
                obj.E = interp1(obj.thermallibs{1}.kT,Evec,E/obj.kT,...
                    'linear','extrap')*E;
                
            end


        end
        
        % sample source energy
        function obj = sample_source_energy(obj,opt)
            
            if strcmp(opt,'const')
                
                % assume fixed source at 2.0 MeV
                obj.E = 2.0;
                
            elseif strcmp(opt,'fission');
                
                % sample random number
                rn1 = obj.rng.rand;
            
                % find bin location
                idx = find(obj.chicdf.*(obj.chicdf > rn1),1,'first');
                obj.E = obj.egrid(idx);
                
            else
                error('FATAL==>Source cant be sampled.')
            end
            
        end
        
        % load thermal scattering library
        function obj = load_thermal_lib(obj,name)
            
            % increment number of isotopes
            obj.n_thermal = obj.n_thermal + 1;
            
            % call constructor of xsdata
            obj.thermallibs{obj.n_thermal} = thermal_lib_class(name);
            
            % append to list
            obj.thermal_list{obj.n_thermal} = name;
            
        end
        
        % sample reaction type
        function obj = sample_reaction(obj,E,iso)
           
            % compute absorption xs
            xs_a = sqrt(obj.E_ref/E)*obj.xs_ref;                                                                                                                                          
            
            % get hydrogen scattering xs
            xs_s = interp1(iso.egrid,iso.xs,E,'linear','extrap');
            
            % compute probability of abs
            p = xs_a/(xs_a + xs_s);
            
            % sample random number
            rn = obj.rng.rand;
                        
            % choose reaction type
            if rn < p
                obj.reaction = 'capture';
            else
                obj.reaction = 'scatter';
            end
            
        end
       
    end
    
    methods (Access = private)
    
        % construct cdf
        function obj = watt_fission_cdf(obj)
            
            % load cdf
            disp('Loading fission chi data...')
            load('./pdfdata/fission_chi','E','chicdf')
            
            % set properties
            obj.egrid = E;
            obj.chicdf = chicdf;
                        
        end
        
    end
    
end


classdef pdf_class
    
    % PDF_CLASS contains information about the distribution function
    %   contains all of the random number sampling and pdfs needed by MC
    
    properties (SetAccess = private, GetAccess = public)
      
        rng    % random number generator object
        E      % new sampled energy
        egrid  % energy grid for chi pdf
        chicdf % cdf for chi
        
    end
    
    methods
        
        % constructor
        function obj = pdf_class(seed)
           
            % initialize random number generator with seed
            obj.rng = RandStream('mcg16807','Seed',seed);
            
            % intialize watt fission spectrum cdf
            obj.egrid = linspace(0,20,10000);
            
            % create cdf
            obj = obj.watt_fission_cdf();
            
        end
        
        % sample energy
        function obj = sample_collision_energy(obj,E,alpha)
           
            obj.E = E - E*(1-alpha)*obj.rng.rand;

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
                idx = find(obj.chicdf.*(obj.chicdf >= rn1),1,'first')-1;
                obj.E = obj.egrid(idx);
                
            else
                error('FATAL==>Source cant be sampled.')
            end
            
        end
       
    end
    
    methods (Access = private)
    
        % pdf for watt fission spectrum
        function chi = watt_fission(obj,E)
            
            chi = 0.453*exp(-1.036*E).*sinh(sqrt(2.29*E));
            
        end
        
        % construct cdf
        function obj = watt_fission_cdf(obj)
           
            % allocate
            obj.chicdf = zeros(length(obj.egrid),1);
            
            % create cdf
            for i = 2:length(obj.egrid)
                
                % numerically integrate
                obj.chicdf(i) = trapz(obj.egrid(1:i),obj.watt_fission(obj.egrid(1:i)));
            
            end
            
        end
        
    end
    
end


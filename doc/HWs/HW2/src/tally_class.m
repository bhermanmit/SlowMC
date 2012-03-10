classdef tally_class
    %TALLY_CLASS defines the object for tallying MC quantites
    %   defines the object for tallying MC quantites
    
    properties(SetAccess = private,GetAccess = public)
        
        nbins   % number of energy grid bins
        egrid   % energy grid
        lgrid   % lethargy grid
        counts  % count
        lcounts % letharg bin counts
        eave    % array of average energy in egrid bins
        leave   % arrage of average energy in lethargy bins
        de      % energy spacing
        le      % lethargy spacing
        flux    % flux accumulation
        lflux   % flux in lethargy bins
        xs_ref = 7;
        E_ref = 0.025e-6;
        
    end
    
    methods
        
        % constructor
        function obj = tally_class(nbins,Eo)
            
            % initalize values
            obj.nbins = nbins;
            obj.egrid = linspace(0,Eo,nbins+1);
            obj.eave = linspace((Eo/(nbins*2)),Eo-(Eo/(nbins*2)),nbins);
            obj.counts = zeros(1,nbins);
            obj.flux = zeros(1,nbins);
            obj.de = obj.egrid(2) - obj.egrid(1);
            
            obj.lgrid = logspace(-10,log10(Eo),nbins+1);
            obj.leave = 0.5*(obj.lgrid(1:nbins) + obj.lgrid(2:nbins+1));
            obj.lcounts = zeros(1,nbins);
            obj.lflux = zeros(1,nbins);
            obj.le = log(max(obj.lgrid)/obj.lgrid(1)) - ...
                     log(max(obj.lgrid)/obj.lgrid(2));
            
        end
        
        % bank tally
        function obj = bank_tally(obj,E,totxs)
           
            % find out bin index
            idx = find(obj.egrid.*(obj.egrid > E),1,'first')-1;
            
            if idx > length(obj.egrid) || idx == 0
                error('FATAL ==> index out of bounds');
            end
            
            % bank a count
            obj.counts(idx) = obj.counts(idx) + 1/totxs;
            
            % find out lethargy bin index
            idx = find(obj.lgrid.*(obj.lgrid >= E),1,'first')-1;
                       
            % bank a count (disregard below the cutoff energy)
            if (idx > 0)
                obj.lcounts(idx) = obj.lcounts(idx) + 1/totxs;
            end
            
        end
        
        % plot tally
        function plot_tally(obj,n)
            
            % get random rgb code
            r = rand(1);
            g = rand(1);
            b = rand(1);
            
            figure(n);
            hold on;
            plot(obj.eave,obj.counts,'Color',[r,g,b]);
            
        end
        
        % plot flux
        function plot_flux(obj,n)
            
            %figure(n)
            %loglog(obj.eave,obj.flux);
            
            figure(n)
            loglog(obj.leave,obj.lflux);
            drawnow;
        end
        
        % clear tally
        function obj = reset(obj)
           
            % bank in flux first
            obj.flux = obj.flux + obj.counts;
            obj.lflux = obj.lflux + obj.lcounts;
            
            % now reset
            obj.counts(:) = 0;
            obj.lcounts(:) = 0;
            
        end
        
    end
    
    methods(Access = private)
       
        % determine the index on the energy grid
        function idx = get_energy_idx(obj,E)
           
            idx = ceil(E/obj.de);
            
        end
        
        % determine the index on the lethargy grid
        function idx = get_lethargy_idx(obj,E)
            
            % this is complex since we are putting on a log energy grid
            l = log(max(obj.lgrid)/E);
            idx = length(obj.lgrid) - floor(l/obj.le) - 1;
            if idx <= 0
                idx = 1;
            end
            
        end
        
    end
    
end


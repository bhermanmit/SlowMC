% Resonance Integral Test

% load U-238 cross sections
load('U_238');

% extract capture cross section
idx = find(strcmp(MT,'abs'),1,'first');
xs_abs = xs(:,idx);

% find max energy 
Etop = E(length(E));

% set energy range
bin = [0.01,0.1,1.0,6.0,10.0,25.0,50.0,100.0,Etop];

% preallocate vectors
RI = zeros(length(bin));


% begin loop over energy range
for i = 1:length(bin) - 1
    
    % energy index bounds
    idxl = find(E>bin(i),1,'first');
    idxr = find(E>bin(i+1),1,'first');
    
    % perform numerical integration
    RI(i) = trapz(E(idxl:idxr),xs_abs(idxl:idxr)./E(idxl:idxr));
    
    % display info
    fprintf('Bin %8.1e - %8.1e:  Resonance Integral: %8.2f\n',bin(i),bin(i+1),RI(i));
    
end

% do whole resonance energy range
idxl = find(E>0.5,1,'first');
idxr = length(E);
RI(length(bin)) = trapz(E(idxl:idxr),xs_abs(idxl:idxr)./E(idxl:idxr));

% display info
fprintf('Bin %8.1e - %8.1e:  Resonance Integral: %8.2f\n',0.5,10000,RI(length(bin)));
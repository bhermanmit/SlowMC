% calculation of cdf for fission energy spectrum

% size
size = 10000;

% energy space
E = linspace(0,20,size);

% create cdf
chicdf = zeros(1,size);

% function for fission
chi = 0.453*exp(-1.036*E).*sinh(sqrt(2.29*E));

% loop around and integrate
for i = 2:size
    
    % integrate
    chicdf(i) = trapz(E(1:i),chi(1:i));
    
end

save('fission_chi','E','chicdf');
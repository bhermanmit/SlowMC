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

% transform to uniform cdf
size = 100000;
chicdf_uni = linspace(0,1,size);
E_uni = interp1(chicdf,E,chicdf_uni);

% write out hdf5 file
delete('fission.h5');
h5create('fission.h5','/vecsize',1);
h5write('fission.h5','/vecsize',size);
h5create('fission.h5','/E',size);
h5write('fission.h5','/E',E_uni);
h5create('fission.h5','/cdf_width',1);
h5write('fission.h5','/cdf_width',chicdf_uni(2)-chicdf_uni(1));
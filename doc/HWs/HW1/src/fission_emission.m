%%  Verification of Fission Emission Spectrum

% brute force the tallying here since we are just verifying the sampling
n_samples = 10000;
seed = 25;
source = 'fission';


% get pdf object
pdf = pdf_class(seed);

% bin structure
nbins = 1000;
egrid = linspace(0,20.0,nbins+1);
de = egrid(2) - egrid(1);
eave = linspace((20.0/(nbins*2)),20.0-(20.0/(nbins*2)),nbins);
counts = zeros(1,nbins);

% begin loop
for j = 1:n_samples

    % sample energy
    pdf = pdf.sample_source_energy(source);
    
    % get index to bin
    idx = find(egrid.*(egrid >= pdf.E),1,'first')-1;
    
    % bank count
    counts(idx) = counts(idx) + 1;
    
    if mod(j,1000) == 0
        fprintf('Samples: %d ...\n',j);
    end
    
end

plot(eave,counts)
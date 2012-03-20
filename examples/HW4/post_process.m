% Post Processing Script for SlowMC code

% get reaction rates
flux    = h5read('output.h5',horzcat('/tally_',num2str(1),'/mean'));
totRR   = h5read('output.h5',horzcat('/tally_',num2str(2),'/mean'));
absRR   = h5read('output.h5',horzcat('/tally_',num2str(3),'/mean'));
scatRR  = h5read('output.h5',horzcat('/tally_',num2str(4),'/mean'));
fissRR  = h5read('output.h5',horzcat('/tally_',num2str(5),'/mean')); 
nfissRR = h5read('output.h5',horzcat('/tally_',num2str(6),'/mean'));

% compute cross sections
totxs   = sum(totRR,2)./sum(flux,2);
absxs   = sum(absRR,2)./sum(flux,2);
scatxs  = sum(scatRR,2)./sum(flux,2);
fissxs  = sum(fissRR,2)./sum(flux,2);
nfissxs = sum(nfissRR,2)./sum(flux,2);

% compute kinf from xs
kinf = sum(sum(nfissRR)) / sum(sum(absRR));

% extract flux finer distribution
flux2 = h5read('output.h5',horzcat('/tally_',num2str(7),'/mean'));

% compute flux ratio mod/fuel
flux_rat = flux2(:,2)./flux2(:,1);

% extract energy and spectrum mean
E_spec = logspace(-11,log10(20.0),1000);
mean_spec = h5read('output.h5',horzcat('/tally_',num2str(8),'/mean'));

% plot spectrum
loglog(E_spec,mean_spec(:,1))
hold on
loglog(E_spec,mean_spec(:,2),'r')
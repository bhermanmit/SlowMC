% Post Processing Script for SlowMC code

% get reaction rates
flux    = h5read('output.h5',horzcat('/tally_',num2str(1),'/mean'));
totRR   = h5read('output.h5',horzcat('/tally_',num2str(2),'/mean'));
absRR   = h5read('output.h5',horzcat('/tally_',num2str(3),'/mean'));
scatRR  = h5read('output.h5',horzcat('/tally_',num2str(4),'/mean'));
fissRR  = h5read('output.h5',horzcat('/tally_',num2str(5),'/mean')); 
nfissRR = h5read('output.h5',horzcat('/tally_',num2str(6),'/mean'));
diffRR  = h5read('output.h5',horzcat('/tally_',num2str(7),'/mean'));
transRR = h5read('output.h5',horzcat('/tally_',num2str(8),'/mean'));

% fast to thermal
fast_to_thermal = sum(flux(2,:))/sum(flux(1,:));

% compute cross sections
totxs   = sum(totRR,2)./sum(flux,2);
absxs   = sum(absRR,2)./sum(flux,2);
scatxs  = sum(scatRR,2)./sum(flux,2);
fissxs  = sum(fissRR,2)./sum(flux,2);
nfissxs = sum(nfissRR,2)./sum(flux,2);
transxs = sum(transRR,2)./sum(flux,2);

% compute estimates of diffusion coefficient
diff_iso = 1./(3*totxs);
diff_trans = 1./(3*transxs);
diffcof = sum(diffRR,2)./sum(flux,2);

% compute kinf from xs
kinf = sum(sum(nfissRR)) / sum(sum(absRR));

% compute effective group 1->2 scattering xs
scat12xs = (nfissxs(2)/kinf - absxs(2))/(1 - nfissxs(1)/(kinf*absxs(1)));

% extract flux finer distribution
flux2 = h5read('output.h5',horzcat('/tally_',num2str(9),'/mean'));

% compute flux ratio mod/fuel
flux_rat = flux2(:,2)./flux2(:,1);

% extract energy and spectrum mean
E_spec = logspace(-11,log10(20.0),5000);
mean_spec = h5read('output.h5',horzcat('/tally_',num2str(10),'/mean'));

% plot spectrum
loglog(E_spec,mean_spec(:,1))
hold on
loglog(E_spec,mean_spec(:,2),'r')

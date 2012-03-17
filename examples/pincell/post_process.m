% Post Processing Script for SlowMC code

% select which tally to read
spectrum = 2;

% extract energy and spectrum mean
E_spec = logspace(-11,log10(20.0),1000);
mean_spec = h5read('output.h5',horzcat('/tally_',num2str(spectrum),'/mean'));

% plot spectrum
loglog(E_spec,mean_spec(:,1))
hold on
loglog(E_spec,mean_spec(:,2),'r')
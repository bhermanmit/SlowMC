% Post Processing Script for SlowMC code

% set tally numbers
user_capt = 1;
user_flux = 2;
spectrum = 3;

% make energy sapce for spectrum
E_spec = logspace(-11,log10(20.0),1000);

% extract flux spectrum
mean_spec = h5read('output.h5',horzcat('/tally_',num2str(spectrum),'/mean'));
std_spec = h5read('output.h5',horzcat('/tally_',num2str(spectrum),'/std'));

% get energy space of user tallies
E_user = h5read('output.h5',horzcat('/tally_',num2str(user_flux),'/E'));

% extract capture tally
mean_capt = h5read('output.h5',horzcat('/tally_',num2str(user_capt),'/mean'));
std_capt = h5read('output.h5',horzcat('/tally_',num2str(user_capt),'/std'));

% extract flux tally
mean_flux = h5read('output.h5',horzcat('/tally_',num2str(user_flux),'/mean'));
std_flux = h5read('output.h5',horzcat('/tally_',num2str(user_flux),'/std'));

% make spectrum plot
loglog(E_spec,mean_spec)

% compute ln factors
LNfactor = zeros(length(E_user) - 1,1);
for i = 1:length(E_user)-1
   LNfactor(i) = log(E_user(i+1)/E_user(i)); 
end

% compute effective resonance integral
RIeff = LNfactor.*(mean_capt./mean_flux);
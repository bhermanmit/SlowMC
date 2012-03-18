% Single Level Breit Wigner xs generator for U-238
tic
% Inputs for user
isoname = 'U_238';  % isotope name
n_res = 14;        % number of resonances to read from file
T = 300;           % temperature of resonances
sig_pot = 11.2934; % potential cross section of isotope
A = 238;           % isotope atomic weight
maxE = 1000;        % max energy in eV

% load isotope file
load('U_238_res.txt');

% constants
k = 8.6173324e-5;

% extract resonance information from loaded dataa
Gg = U_238_res(1:n_res,4);
Gn = U_238_res(1:n_res,3);
E0 = U_238_res(1:n_res,1);

% create energy vector
dE = 0.001;
E = 10.^(log10(1e-5):dE:log10(20e6))';
sizeE = length(E);

% set psi-chi vectors
psi = zeros(sizeE,1);
chi = zeros(sizeE,1);

% initialize xs
xs = zeros(sizeE,2);
xs(:,2) = sig_pot;

% begin loop around resonances
for j = 1:n_res
       
    % psi-chi parameters
    G = Gg(j) + Gn(j);
    r = 2603911/E0(j)*((A+1)/A);
    q = sqrt(r*sig_pot);
    xi = G*sqrt(A/(4*k*T*E0(j)));
    x = 2*(E-E0(j))/G;

    % compute psi-chi functions
    y = ((x+1i)/2*xi);
    psichi = pi*xi/(2*sqrt(pi))*W(y); % compute complex value
    psi = real(psichi);
    chi = 2*imag(psichi);

    % compute xs
    xs(:,1) = (Gn(j)/G)*(Gg(j)/G)*sqrt(E0(j)./E).*(r*psi) + xs(:,1);
    xs(:,2) = (Gn(j)/G)*(Gn(j)/G)*(r*psi + q*chi) + xs(:,2);
    
    % display resonance info
    fprintf('Completed Resonance: %d\n',E0(j));
    
end

% append pseudo resonances (25eV spacing)
Elast = E0(length(E0));
clear E0;
E0 = Elast + 25;
Gg = 0.023;

% begin loop around building pseudo resonances
while E0 < maxE
   
    % compute neutron width
    Gn = 0.050*sqrt(E0/Elast);
    
    % psi-chi parameters
    G = Gg + Gn;
    r = 2603911/E0*((A+1)/A);
    q = sqrt(r*sig_pot);
    xi = G*sqrt(A/(4*k*T*E0));
    x = 2*(E-E0)/G;

    % compute psi-chi functions
    y = ((x+1i)/2*xi);
    psichi = pi*xi/(2*sqrt(pi))*W(y); % compute complex value
    psi = real(psichi);
    chi = imag(psichi);

    % compute xs
    xs(:,1) = (Gn/G)*(Gg/G)*sqrt(E0./E).*(r*psi) + xs(:,1);
    xs(:,2) = (Gn/G)*(Gn/G)*(r*psi + q*chi) + xs(:,2);
    
    % get next E0
    E0 = E0 + 25;
    
    % display resonance info
    fprintf('Completed Resonance: %d\n',E0);
    
end
toc

% change units on E
E = E/1e6;

% % zero out xs over 1 keV
% xs(:,1) = xs(:,1).*(E <= 1e-3);
% 
% % put 0.1 barns after
% xs(:,1) = xs(:,1) + (E > 1e-3)*0.1;
% 
% % get size
% sizeE = length(xs(:,1));
% 
% % get capture
% xs_capt = xs(:,1);
% 
% % set scattering to 0
% xs_scat = sig_pot*ones(sizeE,1);
% xs_fiss = zeros(sizeE,1);
% 
% % filename
% hdfile = horzcat(isoname,'_',num2str(T),'.h5');
% 
% % write out hdf5 file
% delete(hdfile);
% h5create(hdfile,'/vecsize',1);
% h5write(hdfile,'/vecsize',sizeE);
% h5create(hdfile,'/xs_scat',sizeE);
% h5write(hdfile,'/xs_scat',xs_scat);
% h5create(hdfile,'/xs_capt',sizeE);
% h5write(hdfile,'/xs_capt',xs_capt);
% h5create(hdfile,'/xs_fiss',sizeE);
% h5write(hdfile,'/xs_fiss',xs_fiss);
% h5create(hdfile,'/E_width',1);
% h5write(hdfile,'/E_width',dE);
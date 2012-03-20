% H-1 Cross Section

% scattering cross section
run('./H_1_scat.m')

% capture cross section
run('./H_1_capt.m');

% fission cross section
xs_fiss = zeros(length(xs_capt),1);

% get size
sizeE = length(E);

% set up free gas thermal scattering variables
A = 1;
T = 300;
sizeN = 10000;
kTfactor = [0.01,0.05,0.1,0.25,0.5,0.75,1,2,5,10,15,20,25,30,40,50,75,100,125,150,200,300];
[thermalcdf,Erat]  = free_gas(A,T,sizeN,kTfactor);


% write out hdf5 file
hdfile = 'H_1.h5';
delete(hdfile);
h5create(hdfile,'/vecsize',1);
h5write(hdfile,'/vecsize',sizeE);
h5create(hdfile,'/xs_scat',sizeE);
h5write(hdfile,'/xs_scat',xs_scat);
h5create(hdfile,'/xs_capt',sizeE);
h5write(hdfile,'/xs_capt',xs_capt);
h5create(hdfile,'/xs_fiss',sizeE);
h5write(hdfile,'/xs_fiss',xs_fiss);
h5create(hdfile,'/E_width',1);
h5write(hdfile,'/E_width',dE);

h5create(hdfile,'/kTsize',1);
h5write(hdfile,'/kTsize',length(kTfactor));

h5create(hdfile,'/cdfsize',1);
h5write(hdfile,'/cdfsize',length(thermalcdf));

h5create(hdfile,'/kT',length(kTfactor));
h5write(hdfile,'/kT',kTfactor);

h5create(hdfile,'/thermalcdf',length(thermalcdf));
h5write(hdfile,'/thermalcdf',thermalcdf);

h5create(hdfile,'/Erat',[size(Erat,1),size(Erat,2)]);
h5write(hdfile,'/Erat',Erat);

h5create(hdfile,'/cdf_width',1);
h5write(hdfile,'/cdf_width',thermalcdf(2) - thermalcdf(1));
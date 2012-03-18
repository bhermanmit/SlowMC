% U-235 Cross Sections

% Load Capture
run('./U_235_capt.m');

% Load Fission
run('./U_235_fiss.m');

% Make Potential Scattering
xs_scat = 11.5860*ones(length(xs_capt),1);

sizeE = length(E);
hdfile = 'U_235.h5';

% write out hdf5 file
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
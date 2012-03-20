% O-16 Cross Sections

% run scattering
run('./O_16_scat.m');

% set capture and fission to zero
xs_capt = zeros(length(xs_scat),1);
xs_fiss = zeros(length(xs_scat),1);

% size of vectors
sizeE = length(E);

% write out hdf5 file
hdfile = 'O_16.h5';
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
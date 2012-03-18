% U-238 Cross Sections

% run fission
run('./U_238_fiss.m');

% run capture
run('./U_238_capt.m');

% put scattering
xs_scat = 11.2934*ones(length(xs_capt),1);

% SLBW
run('./SLBW.m');

% set U-238 XS on top of U_238 capture
idx = find((E>1e-3),1,'first');
xs_capt(1:idx) = xs(1:idx,1);

% write out HDF5 file
hdfile = horzcat(isoname,'_',num2str(T),'.h5');
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

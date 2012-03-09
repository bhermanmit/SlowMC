% 1/v absorber file

% inputs
xs0 = 2;
E0 = 0.0253*1e-6;

% E range
dE = 0.001;
E = 10.^(log10(1e-11):dE:log10(20.0))';

% make cross sections
xs_capt = sqrt(E0./E)*xs0;
xs_scat = zeros(length(xs_capt),1);

% get size
sizeE = length(xs_capt);

% write out hdf5 file
delete('v_abs.h5');
h5create('v_abs.h5','/vecsize',1);
h5write('v_abs.h5','/vecsize',sizeE);
h5create('v_abs.h5','/xs_scat',sizeE);
h5write('v_abs.h5','/xs_scat',xs_scat);
h5create('v_abs.h5','/xs_capt',sizeE);
h5write('v_abs.h5','/xs_capt',xs_capt);
h5create('v_abs.h5','/E_width',1);
h5write('v_abs.h5','/E_width',dE);
E_endf = [
         1e-11
   1.10304e-11
   1.21669e-11
   1.34205e-11
   1.48034e-11
   1.63287e-11
   1.80111e-11
   1.98669e-11
   2.19139e-11
   2.41719e-11
   2.66625e-11
   2.94097e-11
     3.244e-11
   3.57825e-11
   3.94694e-11
   4.35362e-11
   4.80221e-11
   5.28819e-11
   5.83307e-11
   6.43409e-11
   7.09704e-11
    7.8283e-11
    8.6349e-11
   9.52462e-11
    1.0506e-10
   1.15885e-10
   1.27826e-10
   1.40996e-10
   1.55524e-10
   1.71549e-10
   1.89225e-10
   2.08722e-10
   2.30228e-10
    2.5395e-10
   2.80116e-10
   3.08979e-10
   3.40815e-10
   3.75931e-10
   4.14666e-10
   4.57392e-10
    5.0452e-10
   5.56504e-10
   6.13845e-10
   6.77094e-10
   7.46859e-10
   8.23813e-10
   9.08697e-10
    1.00233e-9
     1.1056e-9
    1.21952e-9
    1.34518e-9
    1.48378e-9
    1.63666e-9
     1.8053e-9
    2.00536e-9
    2.22748e-9
    2.47432e-9
    2.74839e-9
    3.05296e-9
    3.39112e-9
    3.76691e-9
    4.18415e-9
    4.64783e-9
    5.16263e-9
    5.73475e-9
    6.36995e-9
    7.07585e-9
     7.8596e-9
    8.73058e-9
    9.69761e-9
    1.08477e-8
    1.21343e-8
    1.35734e-8
    1.51832e-8
    1.69839e-8
    1.91321e-8
    2.15511e-8
    2.42771e-8
    2.52209e-8
       2.53e-8
    2.86998e-8
    3.25565e-8
     3.7192e-8
    4.24854e-8
    4.88745e-8
    5.66211e-8
    6.60549e-8
    7.76004e-8
    9.24546e-8
    1.11701e-7
    1.36864e-7
    1.71253e-7
    2.20371e-7
    2.93678e-7
     4.1105e-7
    6.12785e-7
     1.0076e-6
    1.96022e-6
    5.33729e-6
    4.33445e-5
    0.00973874
     0.0183791
     0.0246658
     0.0408441
     0.0608465
     0.0811287
      0.100705
      0.151058
           0.2
           0.3
           0.4
           0.5
      0.604114
      0.704768
      0.805643
      0.906349
       1.00705
       1.20846
       1.40988
       1.61129
        1.8127
       2.02822
       2.23104
       2.43386
       2.63668
       2.85953
       3.06379
       3.29093
       3.54594
        3.8808
           4.2
       4.49355
       4.86772
       5.14208
        5.5388
       6.04232
       6.54585
       7.09876
       7.60581
       8.11287
       8.61992
       9.12697
       9.63403
       10.2126
       10.7984
          11.5
       12.2551
       13.0917
            14
       15.1058
       16.1129
       17.1199
        18.127
        19.134
            20
];

xs_endf = [
1159.67
 1104.2
1051.38
1001.08
953.199
907.606
864.195
822.863
783.509
 746.04
710.364
676.397
644.057
613.265
583.948
556.036
 529.46
504.578
480.468
457.513
435.659
414.852
395.043
376.184
358.229
341.136
324.864
309.372
294.625
280.586
267.221
  254.5
 242.39
230.863
219.892
209.449
 199.51
 190.05
181.048
172.481
164.329
156.571
149.191
142.169
135.488
129.134
 123.09
117.342
111.876
106.679
101.738
97.0419
92.5788
88.3381
84.0287
79.9529
76.0955
72.4494
69.0011
65.7442
62.6665
59.7623
57.0207
54.4365
     52
49.7065
47.5472
45.5179
43.6108
41.8219
40.0367
38.3727
 36.824
35.3848
34.0497
32.7392
31.5347
30.4294
30.0957
30.0687
 29.034
28.0958
27.2032
26.4039
25.6551
 24.962
24.3287
 23.757
23.2274
22.7482
22.3237
21.9449
21.6087
 21.316
21.0648
20.8579
20.6926
20.5679
20.4841
20.4368
 19.222
18.2662
17.6757
16.2353
14.8181
13.6515
12.7166
10.8963
9.64353
7.95225
 6.8766
6.12557
5.54899
 5.1153
 4.7652
4.47589
4.23215
3.83739
3.52938
3.27951
3.07085
2.88176
2.72734
2.59131
2.47012
 2.3512
2.25265
2.15322
2.05201
1.93368
1.83279
1.75064
1.65403
1.59141
1.50615
1.41101
 1.3268
  1.245
1.17776
1.11703
1.06191
1.01164
.965605
.917673
.873354
.824917
.778557
.732128
.687115
.639018
.600069
.565181
.533752
.505289
.482757
];

% create new energy space
dE = 0.001;
E = 10.^(log10(1e-11):dE:log10(20.0))';

% create cross sections
xs_scat = interp1(E_endf,xs_endf,E,'linear','extrap');
xs_capt = zeros(length(xs_scat),1);

% get size
sizeE = length(E);

% set up free gas thermal scattering variables
A = 1;
T = 300;
sizeN = 10000;
kTfactor = [0.01,0.05,0.1,0.25,0.5,0.75,1,2,5,10,15,20,25,30,40,50,75,100,125,150,200,300];
[thermalcdf,Erat]  = free_gas(A,T,sizeN,kTfactor);


% write out hdf5 file
delete('H_1.h5');
h5create('H_1.h5','/vecsize',1);
h5write('H_1.h5','/vecsize',sizeE);
h5create('H_1.h5','/xs_scat',sizeE);
h5write('H_1.h5','/xs_scat',xs_scat);
h5create('H_1.h5','/xs_capt',sizeE);
h5write('H_1.h5','/xs_capt',xs_capt);
h5create('H_1.h5','/E_width',1);
h5write('H_1.h5','/E_width',dE);

h5create('H_1.h5','/kTsize',1);
h5write('H_1.h5','/kTsize',length(kTfactor));

h5create('H_1.h5','/cdfsize',1);
h5write('H_1.h5','/cdfsize',length(thermalcdf));

h5create('H_1.h5','/kT',length(kTfactor));
h5write('H_1.h5','/kT',kTfactor);

h5create('H_1.h5','/thermalcdf',length(thermalcdf));
h5write('H_1.h5','/thermalcdf',thermalcdf);

h5create('H_1.h5','/Erat',[size(Erat,1),size(Erat,2)]);
h5write('H_1.h5','/Erat',Erat);

h5create('H_1.h5','/cdf_width',1);
h5write('H_1.h5','/cdf_width',thermalcdf(2) - thermalcdf(1));
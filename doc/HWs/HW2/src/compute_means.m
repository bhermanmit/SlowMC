function compute_means(neut)

% initialize
n_coll_therm = 0;
n_coll = 0;
n_abs = 0;
time = 0;

% sum variables
for i = 1:length(neut)
   
    n_coll_therm = n_coll_therm + neut(i).n_coll_therm;
    n_coll = n_coll + neut(i).n_coll;
    n_abs = n_abs + neut(i).n_abs;
    time = time + neut(i).time;
    
end

% compute means
n_coll_therm = n_coll_therm/length(neut);
n_coll = n_coll/length(neut);
time = time/length(neut);

% display answers
fprintf('Mean number of collisions for neutrons to reach 1.0eV: %d\n',n_coll_therm);
fprintf('Mean number of collisions per neutron: %d\n',n_coll);
fprintf('Mean neutron lifetime [s]: %d\n',time);
fprintf('Number of absorptions per fission neutron: %d\n',n_abs);
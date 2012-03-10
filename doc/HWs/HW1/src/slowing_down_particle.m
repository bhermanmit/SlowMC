% Bryan Herman
% Slowing Down Code
% HW1 22.211
%
%% Slowing down part (by particle)

% set input information
n_histories = 100000;
Eo = 20.0;
A = 1;
seed = 5;
n_bins = 1000;
source = 'fission';

% initalize objects
neut(1:n_histories) = particle_class();
tal = tally_class(n_bins,Eo);
pdf = pdf_class(seed);
mat = material_class(A);

% set materials (may put this in class later)
mat = mat.load_isotope('H_1');

% begin loop around histories
for j = 1:n_histories
    
    % sample source energy
    pdf = pdf.sample_source_energy(source);
    neut(j) = neut(j).set_energy(pdf.E);

    while neut(j).alive == 1
        
        % bank sampled energy (only seed 1 isotope for now)
        tal = tal.bank_tally(neut(j).E,mat.isotopes{1});
        
        % sample new energy
        pdf = pdf.sample_collision_energy(neut(j).E,mat.alpha);
        
        % set particle to that new energy
        neut(j) = neut(j).set_energy(pdf.E);
        
        if(neut(j).E < 1e-8)
            
            % kill particle
            neut(j) = neut(j).kill;
            
            % reset tallies and bank flux
            tal = tal.reset;
            
        end
        
    end
    
    % display calculation progress
    if mod(j,1000) == 0
        fprintf('Histories: %d ...\n',j);
    end
    
end

% plot flux
tal.plot_flux(1);
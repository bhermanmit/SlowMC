% Bryan Herman
% Slowing Down Code
% HW1 22.211
% 
%% Slowing down part (by generation)
% set input information
n_histories = 100000;
n_generations = 15;
Eo = 2.0;
A = 1;
seed = 5;
n_bins = 1000;
source = 'const';

% initalize objects
neut(1:n_histories) = particle_class();
tal = tally_class(n_bins,Eo);
pdf = pdf_class(seed);
mat = material_class(A);

% set materials (may put this in class later)
mat = mat.load_isotope('H_1');

% begin loop
for i = 1:n_generations
    
    % display information
    fprintf('\nExecuting Generation %d\n',i);
    fprintf('=======================\n\n');
    
    for j = 1:n_histories       
        
        % sample source energy if fission spectrum is on
        if i == 1
            pdf = pdf.sample_source_energy(source);
            neut(j) = neut(j).set_energy(pdf.E);
        end
            
        % sample new energy
        pdf = pdf.sample_collision_energy(neut(j).E,mat.alpha);
        
        % set particle to that new energy
        neut(j) = neut(j).set_energy(pdf.E);
            
        % bank sampled energy (only seed 1 isotope for now)
        tal = tal.bank_tally(neut(j).E,mat.isotopes{1});
        
        if mod(j,1000) == 0
            fprintf('Histories: %d ...\n',j);
        end
    end
    
    % plot that generation
    tal.plot_tally(1);
    
    % reset tallies
    tal = tal.reset;
    
end

% plot flux
tal.plot_flux(2);


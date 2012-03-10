% Bryan Herman
% Slowing Down Code
% HW1 22.211
%
%
% clear close everything
clear
close all
clc
profile on
diary('results.txt');
% set input information
n_histories = 1000000;
Eo = 20.0;
A = 1;
seed = 5;
n_bins = 1000;
source = 'fission';

% print header
print_run_info();

% initalize objects
neut(1:n_histories) = particle_class();
tal = tally_class(n_bins,Eo);
pdf = pdf_class(seed);
mat = material_class(A);

% set materials (may put this in class later)
mat = mat.load_isotope('H_1',1);

% load thermal scattering libraries
pdf = pdf.load_thermal_lib('H_1');

% begin loop around histories
fprintf('\nBeginning Histories\n=====================\n')
for j = 1:n_histories
    
    % sample source energy
    pdf = pdf.sample_source_energy(source);
    neut(j) = neut(j).set_energy(pdf.E);

    % perform scattering events until cutoff or absorbed
    while neut(j).alive == 1
        
        % calculate total macro xs
        mat = mat.compute_macro_totxs(neut(j).E);
        
        % bank neutron time
        neut(j) = neut(j).add_time(mat.totxs);
        
        % bank sampled energy (only 1 isotope for now)
        tal = tal.bank_tally(neut(j).E,mat.totxs);
        
        % sample reaction type
        pdf = pdf.sample_reaction(neut(j).E,mat.isotopes{1});
        
        % bank neutron collision information
        neut(j) = neut(j).collision(pdf.reaction);
        
        % check reaction type
        if strcmp(pdf.reaction,'capture')
            
            % kill neutron
            neut(j) = neut(j).kill;
            
        elseif strcmp(pdf.reaction,'scatter')
        
            % sample new energy
            pdf = pdf.sample_collision_energy(neut(j).E,mat.alpha);
        
            % set particle to that new energy
            neut(j) = neut(j).set_energy(pdf.E);
            
        else
            error('FATAL==> could not sample reaction type');
        end
        
        % check if neutron energy is too low (energy cutoff)
        if(neut(j).E < 1e-10)
            
            % kill particle
            neut(j) = neut(j).kill;
            
        end
        
        % check if neutron is dead
        if (neut(j).alive == 0)
            
            % reset tallies and bank flux
            tal = tal.reset;
            
        end
        
    end

    % display calculation progress
    if mod(j,1000) == 0
        fprintf('Histories: %d ...\n',j);
        
        % plot flux
        tal.plot_flux(1);
    end
    
end

% plot flux
tal.plot_flux(1);

profile viewer
diary off

% compute and print averages
compute_means(neut);
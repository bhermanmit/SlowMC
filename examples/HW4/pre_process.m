% Pre-processing script for material number densities

% Avagadros Number
Na = 0.6022;

% Mass Densities
rho_fuel = 10.2;  % UO2
rho_clad = 6.549; % Zr-90
rho_cool = 0.9966; % H2O

% Enrichment (w% U235/w% U)
% enr = 0.03035;
enr = 0.016;

% Molar Masses
M25 = 235.0439231;
M28 = 238.0507826;
M16 = 15.9949146;
M90 = 89.9047037;
MH1 = 1.0078250;

% Calculate Molar Mass of U
MU = (enr/M25 + (1-enr)/M28)^-1;

% Weight percent of Uranium
wU = MU/(MU + 2*M16);

% Number Density of Fuel 
N25 = rho_fuel*(wU*enr)*Na/M25;
N28 = rho_fuel*(wU*(1-enr))*Na/M28;
N16 = rho_fuel*(1-wU)*Na/M16;

% Number Density of Coolant
Mw = 2*MH1 + M16;
NH = rho_cool*2*Na/Mw;
NO = rho_cool*Na/Mw;

% Number Density of Clad
NZr = rho_clad*Na/M90;
% free gas thermal scattering kernel generation

% inputs
isoname = 'H_1';
savelib = 'no';
A = 1;
T = 300;
size = 10000;
kTfactor = [0.01,0.05,0.1,0.25,0.5,0.75,1,2,5,10,15,20,25,30,40,50,75,100,125,150,200,300];

% constants
k = 8.6173324e-5;
eta = (A+1)/(2*sqrt(A));
rho = (A-1)/(2*sqrt(A));
alpha = ((A-1)/(A+1))^2;

% create output space
p = zeros(size,length(kTfactor));
thermalcdf = zeros(size,length(kTfactor));

% begin loop around temperatures
for i = 1:length(kTfactor)
    
    % create energy space
    kT = k*T;
    E = kT*kTfactor(i);
    Ep = linspace(0,5*E,size);
    
    % compute probability
    p(:,i) =((eta^2/2)*(erf(eta*sqrt(Ep/kT) - rho*sqrt(E/kT)) - ...
        erf(eta*sqrt(Ep/kT) + rho*sqrt(E/kT)) + ...
        exp((E-Ep)/kT).*(erf(eta*sqrt(E/kT) - rho*sqrt(Ep/kT)) + ...
        erf(eta*sqrt(E/kT) + rho*sqrt(Ep/kT))))).*(Ep >= E) + ...
        ((eta^2/2)*(erf(eta*sqrt(Ep/kT) - rho*sqrt(E/kT)) + ...
        erf(eta*sqrt(Ep/kT) + rho*sqrt(E/kT)) + ...
        exp((E-Ep)/kT).*(erf(eta*sqrt(E/kT) - rho*sqrt(Ep/kT)) - ...
        erf(eta*sqrt(E/kT) + rho*sqrt(Ep/kT))))).*(Ep < E);
    
    % multiply by alpha
    p(:,i) = p(:,i)*(1-alpha);    
    
    % relative energy ratio
    Erat = Ep/E;
    
    % loop around energy ratio and integerate
    for j = 2:length(Erat)
        
        % perform integral
        thermalcdf(j,i) = trapz(Erat(1:j),p(1:j,i));
        
    end
    
    % renormalize cdf
    thermalcdf(:,i) = thermalcdf(:,i)/thermalcdf(length(thermalcdf),i);
    
    % plot pdfs and cdfs
    rn1 = rand(1,1);
    rn2 = rand(1,1);
    rn3 = rand(1,1);
    figure(3)
    xlabel('Scattered Energy Ratio (out/in)')
    ylabel('Probability')
    hold on
    pdfplot = plot(Erat,p(:,i),'Color',[rn1,rn2,rn3],'LineWidth',2.0);
    if i == 1
        legend(horzcat('E = ',num2str(kTfactor(i)),' kT'))
    else
        [LEGH,OBJH,OUTH,OUTM] = legend;
        legend([OUTH;pdfplot],OUTM{:},horzcat('E = ',num2str(kTfactor(i)),' kT'));
    end
    
    figure(2)
    xlabel('Scattered Energy Ratio (out/in)')
    ylabel('Cumulative Probability')
    hold on
    cdfplot = plot(Erat,thermalcdf(:,i),'Color',[rn1,rn2,rn3],'LineWidth',2.0);
    if i == 1
        legend(horzcat('E = ',num2str(kTfactor(i)),' kT'))
    else
        [LEGH,OBJH,OUTH,OUTM] = legend;
        legend([OUTH;cdfplot],OUTM{:},horzcat('E = ',num2str(kTfactor(i)),' kT'));
    end
    
end

% save output
if strcmp(savelib,'yes')
    kT = kTfactor;
    save(horzcat(isoname,'_therm'),'Erat','kT','thermalcdf');
end
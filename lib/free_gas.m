function [thermalcdf,Erat]  = free_gas(A,T,sizeN,kTfactor)
% free gas thermal scattering kernel generation

% constants
k = 8.6173324e-5;
eta = (A+1)/(2*sqrt(A));
rho = (A-1)/(2*sqrt(A));
alpha = ((A-1)/(A+1))^2;

% create output space
p = zeros(sizeN,length(kTfactor));
thermalcdf = zeros(sizeN,length(kTfactor));

% begin loop around temperatures
for i = 1:length(kTfactor)
    
    % create energy space
    kT = k*T;
    E = kT*kTfactor(i);
    Ep = linspace(0,5*E,sizeN);
    
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

%% Re-adjust cdf for sampling optimization
Enew = zeros(sizeN,length(kTfactor));
cdfnew = linspace(0,1,sizeN)';
thermalcdftemp = thermalcdf;

% correct thermal cdf for precision unique issues
for i = 1:size(thermalcdftemp,2)
    for j = 1:size(thermalcdftemp,1)-1
        if abs(thermalcdftemp(j,i) - thermalcdftemp(j+1,i)) < 1e-8
            thermalcdf(j+1,i) = 1.01*thermalcdf(j,i);
        end
    end
end

% create new energy vector
for i = 1:size(thermalcdf,2)
    Enew(:,i) = interp1(thermalcdf(:,i),Erat,cdfnew,'linear');
end

% bank to saved variables
thermalcdf = cdfnew;
Erat = Enew;
width = thermalcdf(2) - thermalcdf(1);
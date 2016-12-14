function scatter = waveguide1(energy,N0,theta,thetaPrime,dr, yVal, X, File1, File2)
    %Takes the source and models the transport of neutrons in a waveguide to
    %calculate where and how many neutrons scattered in the reflecting material
    %Constants
    R = 56.35; %cm
    L = 150; %cm 
    sigma = 0; %scattering xsection for energy
    scatterHold = zeros(length(dr),1); %place holder for scattering events
    rho = 1.85*10^3; %kg/m3
    M = 9.012*10^(-3); %kg/mol
    sigmaInelastic = 0;
    %sPoint = source/(2*pi*R*h);
    %Equations of waveguide
    yWG = (R+L)/cos(deg2rad(30))+30; %width of waveguide is 30 cm
    %calculate atomic density of Beryllium
    N = (rho*6.022*10^(23)/M)/1000000; %N/cm3
    %Define the scattering cross section of Beryllium
    [EnergyBE, CrossSectionBE] = read_file(File1);
    sigmaBe = CrossSectionBE*10^(-28);
    energyBe = EnergyBE; %eV
    [energyBeInelastic, sigmaBeInelastic] = read_file(File2);
    sigmaBeInelastic = sigmaBeInelastic*10^(-28);
    for i = 1:length(sigmaBe)-1
        if energyBe(length(energyBe)) ~= energy
            if energyBe(i) == energy
                sigma = sigmaBe(i);
            elseif (energyBe(i) < energy) && (energyBe(i+1) > energy)
                sigma = sigmaBe(i) + (energy-energyBe(i))*(sigmaBe(i+1)-sigmaBe(i))/(energyBe(i+1)-energyBe(i));
            end
        else
            sigma = sigmaBe(length(sigmaBe));
        end
    end
    for i = 1:length(sigmaBeInelastic)-1
        if energyBeInelastic(length(energyBeInelastic)) ~= energy
            if energyBeInelastic(i) == energy
                sigmaInelastic = sigmaBeInelastic(i);
            elseif (energyBeInelastic(i) < energy) && (energyBeInelastic(i+1) > energy)
                sigmaInelastic = sigmaBeInelastic(i) + (energy-energyBeInelastic(i))*(sigmaBeInelastic(i+1)-sigmaBeInelastic(i))/(energyBeInelastic(i+1)-energyBeInelastic(i));
            end
        else
            sigmaInelastic = sigmaBeInelastic(length(sigmaBeInelastic));
        end
    end
    %Macroscopic cross section is
    sigMac = sigma*N;
    sigMacInelastic = sigmaInelastic*N;
    %For N0 = incident neutron flux
    %Figure out the scattering position of neutrons
    r = sum(dr); 
    for i = 1:length(dr)
        %determine if scattering events allow particles to reach target
        if yVal*X(i,1) - tan(thetaPrime)*(150-X(i,1)) < (L+R)*tan(30)
            N_scatter = N0*sigMac*exp(-sigma*N*r)*dr(i);
            N_attenuation = N0*sigMacInelastic*exp(-sigmaInelastic*N*r)*dr(i);
            %N_scatter = round(N_scatter);
            N0 = N0 - N_scatter - N_attenuation;
            scatterHold(i) = N_scatter;
            N_scatter = 0;
            N_attenuation = 0;
        end
    end
    scatter = scatterHold;  
end
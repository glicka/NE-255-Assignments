function scatter = waveguide3(energy,N0,theta,thetaPrime,dr, yVal, X, File1, File2)
    %Takes the source and models the transport of neutrons in a waveguide to
    %calculate where and how many neutrons scattered in the reflecting material
    %Constants
    R = 0; %20; %cm
    L = 200; %cm 
    sigma = 0; %scattering xsection for energy
    scatterHold = 0; %place holder for scattering events
    rho = 9.73; %g/cm3
    M = 208.98; %g/mol
    %sPoint = source/(2*pi*R*h);
    %Equations of waveguide
    yWG = -(R+L)/cos(deg2rad(58.66))+30; %width of waveguide is 30 cm
    %calculate atomic density of Bismuth
    N = rho*6.022*10^(23)/M; %N/cm3
    sigmaInelastic = 0;
    %Define the scattering cross section of Bismuth
    [EnergyBi, CrossSectionBi] = read_file(File1);
    sigmaBi = CrossSectionBi*10^(-28); %cm2
    energyBi = EnergyBi; %eV
    [energyBiInelastic, sigmaBiInelastic] = read_file(File2);
    sigmaBiInelastic = sigmaBiInelastic*10^(-28);
    for i = 1:length(sigmaBi)-1
        if energyBi(length(energyBi)) ~= energy
            if energyBi(i) == energy
                sigma = sigmaBi(i);
            elseif (energyBi(i) < energy) && (energyBi(i+1) > energy)
                sigma = sigmaBi(i) + (energy-energyBi(i))*(sigmaBi(i+1)-sigmaBi(i))/(energyBi(i+1)-energyBi(i));
            end
        else
            sigma = sigmaBi(length(sigmaBi));
        end
    end
    for i = 1:length(sigmaBiInelastic)-1
        if energyBiInelastic(length(energyBiInelastic)) ~= energy
            if energyBiInelastic(i) == energy
                sigmaInelastic = sigmaBiInelastic(i);
            elseif (energyBiInelastic(i) < energy) && (energyBiInelastic(i+1) > energy)
                sigmaInelastic = sigmaBiInelastic(i) + (energy-energyBiInelastic(i))*(sigmaBiInelastic(i+1)-sigmaBiInelastic(i))/(energyBiInelastic(i+1)-energyBiInelastic(i));
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
        
        N_scatter = N0*sigMac*exp(-sigma*N*r)*dr(i);
        N_attenuation = N0*sigMacInelastic*exp(-sigmaInelastic*N*r)*dr(i);
        %N_scatter = round(N_scatter);
        N0 = N0 - N_scatter - N_attenuation;
        scatterHold = scatterHold+N_scatter;

    end
    scatter = scatterHold;
end

function flux = diff_neutronflux(s, energy, r)
    rho = 5.9; %g/cm3
    N = 6.022*10^23*rho*(1/92.232);
    sigma1 = 0;
    sigma2 = 0;
    sigma3 = 0;
    sigma4 = 0;
    [energyZr, sigmaZr] = read_file('sigZr.txt');
    [energyH, sigmaH] = read_file('sigH.txt');
    [energyZr_I, sigmaA_Zr] = read_file('SigmaA_ZR.txt');
    [energyH_I, sigmaA_H] = read_file('SigmaA_H.txt');
    sigmaZr = sigmaZr * 10^(-28);
    sigmaH = sigmaH * 10^(-28);
    sigmaA_Zr = sigmaA_Zr * 10^(-28);
    sigmaA_H = sigmaA_H * 10^(-28);

    %We need to import the sigma values for Zr and H
    %Calculate sigma for ZrH
    for n = 1:length(energy)
        for i = 1:length(sigmaZr)-1
            if energyZr(length(energyZr)) ~= energy(n)
                if energyZr(i) == energy(n)
                    sigma1 = sigmaZr(i);
                elseif (energyZr(i) < energy(n)) && (energyZr(i+1) > energy(n))
                    sigma1 = sigmaZr(i) + (energy(n)-energyZr(i))*(sigmaZr(i+1)-sigmaZr(i))/(energyZr(i+1)-energyZr(i));
                end
            else
                sigma1 = sigmaZr(length(sigmaZr));
            end
        end
        for i = 1:length(sigmaH)-1
            if energyH(length(energyH)) ~= energy(n)
                if energyH(i) == energy(n)
                    sigma2 = sigmaH(i);
                elseif (energyH(i) < energy(n)) && (energyH(i+1) > energy(n))
                    sigma2 = sigmaH(i) + (energy(n)-energyH(i))*(sigmaH(i+1)-sigmaH(i))/(energyH(i+1)-energyH(i));
                end
            else
                sigma2 = sigmaH(length(sigmaH));
            end
        end
        for i = 1:length(sigmaA_Zr)-1
            if energyZr_I(length(energyZr_I)) ~= energy(n)
                if energyZr_I(i) == energy(n)
                    sigma3 = sigmaA_Zr(i);
                elseif (energyZr_I(i) < energy(n)) && (energyZr_I(i+1) > energy(n))
                    sigma3 = sigmaA_Zr(i) + (energy(n)-energyZr_I(i))*(sigmaA_Zr(i+1)-sigmaA_Zr(i))/(energyZr_I(i+1)-energyZr_I(i));
                end
            else
                sigma3 = sigmaA_Zr(length(sigmaA_Zr));
            end
        end
        for i = 1:length(sigmaA_H)-1
            if energyH_I(length(energyH_I)) ~= energy(n)
                if energyH_I(i) == energy(n)
                    sigma4 = sigmaA_H(i);
                elseif (energyH_I(i) < energy(n)) && (energyH_I(i+1) > energy(n))
                    sigma4 = sigmaA_H(i) + (energy(n)-energyH_I(i))*(sigmaA_H(i+1)-sigmaA_H(i))/(energyH_I(i+1)-energyH_I(i));
                end
            else
                sigma4 = sigmaA_H(length(sigmaA_H));
            end
        end
        mac_cross_section = N*(.978*sigma1 + .022*sigma2);
        D = 1/(3*mac_cross_section); %ZrH2 constant
        absorption_cross_section = N*(.978*sigma3 + .022*sigma4);
        L = sqrt(D/absorption_cross_section);
        flux = (s/D) * exp(-r/L);
    end
end
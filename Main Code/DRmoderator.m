function DRmoderator
    doseRate = [2.63;
        2.7567;
        5.0705;
        8.518;
        11.556;
        17.498;
        27.168;
        32.228;
        33.599;
        32.515;
        31.788;
        30.974;
        30.057;
        27.168];
    energy = [2.10E+06;
        2.42E+06;
        2.74E+06;
        3.65E+06;
        4.55E+06;
        5.62E+06;
        7.31E+06;
        9.09E+06;
        1.14E+07;
        1.33E+07;
        1.42E+07;
        1.51E+07;
        1.65E+07;
        1.87E+07];
    percent = [0.008960018;
        0.009391666;
        0.017274438;
        0.029019556;
        0.039369569;
        0.059613077;
        0.092557325;
        0.109795991;
        0.114466784;
        0.110773757;
        0.108296977;
        0.1055238;
        0.102399718;
        0.092557325];
    sigmaTot = zeros(length(doseRate),1);
    x = 0;
    newDoseRate = 0;
    rho = 5.9; %g/cm3
    N = 6.022*10^23*rho*(1/92.232);
    sigma1 = 0;
    sigma2 = 0;
    [energyZr, sigmaZr] = read_file('sigZr.txt');
    [energyH, sigmaH] = read_file('sigH.txt');
    sigmaZr = sigmaZr * 10^(-28);
    sigmaH = sigmaH * 10^(-28);
    trueX = 0;
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
        sigma1;
        sigma2;
        sigmaTot(n) = 0.5*(sigma1 + sigma2);
        %lambda = 10^(-10)*(4.1357*10^(15)*3*10^8)/energy(n);
        %sigmaTot(n) = (1/N)*(1.2865+1.5188*lambda-0.13603*lambda^2+0.0052822*lambda^3);

    end
    for i = 1:length(doseRate)
        x = -log(2/doseRate(i))/(sigmaTot(i)*N);
        trueX = trueX + 1.9*percent(i)*x;
    end
    for i = 1:length(doseRate)
        newDoseRate = newDoseRate + doseRate(i)*exp(-sigmaTot(i)*N*trueX);
    end
    trueX = trueX*10^(-2);
    display(trueX)
    display(newDoseRate)
end
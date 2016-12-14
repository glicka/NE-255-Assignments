function source2 = n2n(N0,R0,energy,scatterMatrix1,scatterMatrix2, File1, File2)
    %Constants
    n2nMatrix = zeros(length(energy),5);
    [rho, M] = read_file(File2); %kg/m3    %kg/mol
    theta = scatterMatrix1(:,2);
    phi = scatterMatrix2(:,2);
    N1 = scatterMatrix1(:,1);
    N2 = scatterMatrix2(:,1);
    sigma = 0;
    %calculate atomic density of Beryllium
    N = rho*6.022*10^(23)/(M*1000000); %N/cm3
    [energy_n2nBe, sigma_n2nBe] = read_file(File1);
    sigma_n2nBe = sigma_n2nBe * 10^(-28); %cm2 
    for k = 1:length(phi)
        %Calculate the (n,2n) sigma for the given energy
        if energy > 1.7*10^(6)
            for i = 1:length(sigma_n2nBe)-1
                if energy_n2nBe(length(energy_n2nBe)) ~= energy
                    if energy_n2nBe(i) == energy
                        sigma = sigma_n2nBe(i);
                    elseif (energy_n2nBe(i) < energy) && (energy_n2nBe(i+1) > energy)
                        sigma = sigma_n2nBe(i) + (energy-energy_n2nBe(i))*(sigma_n2nBe(i+1)-sigma_n2nBe(i))/(energy_n2nBe(i+1)-energy_n2nBe(i));
                    end
                else
                    sigma = sigma_n2nBe(length(sigma_n2nBe));
                end
            end
            %Calculate the thickness of material particles pass through
            r1 = R0/cos(deg2rad(theta(k)));
            r2 = R0/cos(deg2rad(phi(k)));
            %Calculate the number of neutrons that are produced by (n,2n)
            %reactions
            N_n2n = 2*N1(k)*(1-exp(-sigma*N*r1))+2*N2(k)*(1-exp(-sigma*N*r2));
            %Calculate the energy of resulting neutrons
            energy_n2n = 0.43*energy;
            N_pass = N1(k)*exp(-sigma*N*r1)+N2(k)*exp(-sigma*N*r2);
            %Form the matrix with neutrons generated in the target and
            %neutrons that passed through without interacting
            n2nMatrix(k,:) = [N_n2n energy_n2n N_pass 0 0];
        end
    end
    n2nMatrix(1,4) = 2*N0*(1-exp(-sigma*N*R0));
    n2nMatrix(1,5) = N0*exp(-sigma*N*R0);
    source2 = n2nMatrix;
end



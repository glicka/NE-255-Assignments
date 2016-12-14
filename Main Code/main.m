function main(EnergyFlux, FluxSource, File1, File2, File3, File4, File5, File6)
    %Given an input source, models the neutron transport in a neutron treatment
    %waveguide to deliver a beam for cancer therapy.
    %Constants
    phi0 = 47.9; %degrees
    theta0 = 30; %degrees
    thetaNew0 = 11.31; %degrees
    phiNew0 = 11.31; %degrees
    R = 56.35; %cm
    L = 150; %cm
    L1 = 200; %cm
    x = 0;
    doseRate = 0;
    temp = 0;
    temp1 = 0;
    doseRate1 = 0;
    doseRate3 = 0;
    SLD = 9.6301*10^(-06); %scattering length density for Be
    dr = zeros(93,1); %differential path length
    scatterCount1 = 1;
    scatterCount2 = 1;
    fluxEnd = 0;
    scatterMatrix1 = zeros(36,2); %[N_scatter, thetaPrime]
    scatterMatrix2 = zeros(36,2); %[N_scatter, phiPrime]
    scatter3 = 0;
    scatter4 = 0;
    sIncident = 0;

    %Create meshgrid
    x = 1:83; %2.5 cm
    [X,Y] = meshgrid(x);
    y = 1:92; %2.5 cm
    [A,B] = meshgrid(y);

    %Start with energy discretization
    %These values taken from Serpent
    energyFlux = energyOut(EnergyFlux); %MeV
    energyFlux = energyFlux*10^(6); %MeV to eV
    fluxSource = fluxOut(FluxSource);

    %Convert flux into point source
    area = 125*25*5; %area of a hexagon
    sPoint = -1*fluxSource/area;
    doseRateTest = zeros(length(sPoint),1);
    doseRateTest1 = zeros(length(sPoint),1);
    doseRateTest3 = zeros(length(sPoint),1);
    for i = 1:length(sPoint)

        if energyFlux(i) > 1.7*10^6
            %Calc particle wavelength for the energy
            lambda = (4.1357*10^(15)*3*10^8)/energyFlux(i); %hc/E = eV*s*m/s/eV = m
            %Calc neutron refraction index
            nSnell = 1 - lambda^2*SLD/(2*pi); %MATLAB rounds to 1 for most cases
            %Then we visit angular discretization
            thetas = linspace(30,43.3,6);
            sweepTheta = linspace(0,47.7,6);
            for theta = thetas %degrees
                thetaRot = theta-theta0; %Changes reference frame
                thetaPrime = acos((1/nSnell)*cos(deg2rad(thetaRot)));
                thetaPrime = rad2deg(thetaPrime);
                thetaPrime = thetaPrime + theta0; %outgoing total angle
                x1 = solving_lines1(theta);
                %x1 = x1*20;
                y = tan(deg2rad(theta))*x1;
                for indexTheta = sweepTheta
                    dr = zeros(83,1);
                    %Calculate the number of neutrons scattered at
                    %(thetaPrime)
                    y2 = solving_lines2(indexTheta);
                    %y2 = y2*20;
                    z = tan(deg2rad(indexTheta))*y2;
                    for count = 1:83
                        if count > 22 && count < x1/2.5
                            dr(count) = sqrt((X(1,count)-X(1,count-1))^2 + (y*X(1,count)-y*X(1,count-1))^2 + (z*X(1,count)-z*X(1,count-1))^2);
                        end
                    end
                    dr = (L/(cos(deg2rad(theta))*cos(deg2rad(indexTheta))))*dr/sum(dr);
                    scatter1 = waveguide1(energyFlux(i),sPoint(i),indexTheta,thetaPrime,dr, y, X, File1, File2);
                    diffuseFlux = diffusionZrH2(scatter1,dr,X,y,thetaPrime, energyFlux(i));
                    dFlux = sum(diffuseFlux);
                    scatterMatrix1(scatterCount1,:) = [dFlux thetaPrime];
                    scatterCount1 = scatterCount1+1;
                end
            end
            phis = linspace(42,47.7,6);
            sweepPhi = linspace(0,43.3,6);
            for phi = phis %degrees
                phiRot = phi-phi0; %Changes reference frame
                %Determines the exit angle from interaction
                phiPrime = acos((1/nSnell)*cos(deg2rad(phiRot)));
                phiPrime = rad2deg(phiPrime);
                phiPrime = phiPrime + phi0; %outgoing total angle
                %Calculate the number of neutrons scattered at
                %(thetaPrime,phiPrime)
                x2 = solving_lines2(phi);
                y = tan(deg2rad(phi))*x2;
                for indexPhi = sweepPhi
                    dr = zeros(83,1);
                    y2 = solving_lines1(indexPhi);
                    z = tan(deg2rad(indexPhi))*y2;
                    for count = 1:83
                        if count > 22 && count < x2/2.5
                            dr(count) = sqrt((X(1,count)-X(1,count-1))^2 + (y*X(1,count)-y*X(1,count-1))^2 + (z*X(1,count)-z*X(1,count-1))^2);
                        end
                    end
                    dr = (L/(cos(deg2rad(phi))*cos(deg2rad(indexPhi))))*dr/sum(dr);
                    scatter2 = waveguide2(energyFlux(i),sPoint(i),indexPhi,phiPrime,dr, y, X, File1, File2);
                    diffuseFlux = diffusionZrH2(scatter2,dr,X,y2,phiPrime, energyFlux(i));
                    dFlux = sum(diffuseFlux);
                    scatterMatrix2(scatterCount2,:) = [dFlux phiPrime];
                    scatterCount2 = scatterCount2+1;
                end
            end
            %Calculate the neutrons leaving the Beryllium target
            R0=30; %cm (will change if time permits for optimization of thickness)
            % if doseRate < 2 && 0.1*E_total > 2
            % Calculate area of target
            areaBe = (150+82)*tan(30)*(150+82)*tan(47);
            thetaIncident = linspace(0,30,6);
            phiIncident = linspace(0,47.7,6);
            for pleson = 1:6
                for pleson1 = 1:6
                    sIncident = sIncident + diff_neutronflux(sPoint(i),energyFlux(i),L/(cos(deg2rad(thetaIncident(pleson)))*cos(deg2rad(phiIncident(pleson1)))));
                end
            end
           
            sMulti = n2n(sIncident,R0,energyFlux(i),scatterMatrix1,scatterMatrix2, File5, File6);
            sScatterN2N = sMulti(:,1:2)/areaBe; %N/cm2s
            sSourceN2N = sMulti(1,4:5)/areaBe; %N/cm2s
            sPass = sMulti(:,3)/areaBe; %N/cm2s
            sPass = sum(sPass)+sMulti(1,5)/areaBe;
            %Calculate the neutron flux at the end of the waveguide
            %Create a new point source at center of target
            sourceNew = sum(sScatterN2N(:,1)) + sSourceN2N(1,1); %neutrons/s
            
            sourceNew = sourceNew/areaBe; %neutrons/cm2s
            thetasNew = linspace(11.31,90,6);
            thetaSweepNew = linspace(0,49.49,6);
            for thetaNew = thetasNew
                thetaNewRot = thetaNew-thetaNew0;
                thetaNewPrime = acos((1/nSnell)*cos(deg2rad(thetaNewRot)));
                thetaNewPrime = rad2deg(thetaNewPrime);
                thetaNewPrime = thetaNewPrime + thetaNew0; %outgoing total angle
                x3 = solving_lines3(thetaNew);
                y = -tan(deg2rad(thetaNew))*x3;
                for thetaIndexNew = thetaSweepNew
                    dr = zeros(92,1);
                    y2 = solving_lines1(thetaIndexNew);
                    z = -tan(deg2rad(thetaIndexNew))*y2;
                
                for count = 1:92
                    if count > 12 && count < x3/2.5
                        
                        dr(count) = sqrt((A(1,count)-A(1,count-1))^2 + (y*A(1,count)-y*A(1,count-1))^2 + (z*A(1,count)-z*A(1,count-1))^2);
                    end
                end
                if sum(dr)>0
                    dr = (L1/(cos(deg2rad(thetaNew))*cos(deg2rad(thetaIndexNew))))*dr/sum(dr);
                end
                scatter3 = waveguide3(sScatterN2N(1,2),sourceNew,thetaIndexNew,thetaNewPrime,dr, y, A, File3, File4);
                scatter4 = waveguide3(energyFlux(i),sPass,thetaIndexNew,thetaNewPrime,dr, y, A, File3, File4);
                if sScatterN2N(1,2)>1*10^6
                    temp = scatter3/nDoseRate(sScatterN2N(1,2));
                end
                doseRate = doseRate + temp + scatter4/nDoseRate(energyFlux(i));       
                doseRateTest(i) = doseRateTest(i) + temp + scatter4/nDoseRate(energyFlux(i));
                temp = 0;
                end
            end
            phisNew = linspace(11.31,90,6);
            phiSweepNew = linspace(0,58.66,6);
            for phiNew = phisNew
                phiNewRot = phiNew-phiNew0;
                phiNewPrime = acos((1/nSnell)*cos(deg2rad(phiNewRot)));
                phiNewPrime = rad2deg(phiNewPrime);
                phiNewPrime = phiNewPrime + phiNew0; %outgoing total angle
                x4 = solving_lines4(phiNew);
                y = -tan(deg2rad(phiNew))*x4;
                for phiIndexNew = phiSweepNew
                    dr = zeros(92,1);
                    y2 = solving_lines1(phiIndexNew);
                    z = -tan(deg2rad(phiIndexNew))*y2;
                for count = 1:92
                    if count > 12 && count < x4/2.5
                        dr(count) = sqrt((A(1,count)-A(1,count-1))^2 + (y*A(1,count)-y*A(1,count-1))^2 + (z*A(1,count)-z*A(1,count-1))^2);
                    end
                end
                if sum(dr)>0
                    dr = (L1/(cos(deg2rad(phiNew))*cos(deg2rad(phiIndexNew))))*dr/sum(dr);
                end
                scatter5 = waveguide4(sScatterN2N(1,2),sourceNew,phiIndexNew,phiNewPrime,dr, y, A, File3, File4);
                scatter6 = waveguide4(energyFlux(i),sPass,phiIndexNew,phiNewPrime,dr, y, A, File3, File4);
                if sScatterN2N(1,2)>1*10^6
                    temp1 = scatter5/nDoseRate(sScatterN2N(1,2));
                end
                doseRate1 = doseRate1 + temp1 + scatter6/nDoseRate(energyFlux(i));
                doseRateTest1(i) = doseRateTest1(i) + temp1 + scatter6/nDoseRate(energyFlux(i));
                temp1 = 0;
                end
            end
        
        
        %neutrons that directly travel to waveguide exit
        portsTheta = linspace(0,5.7,6);
        portsPhi = linspace(0,5.7,6);
        for portTheta = portsTheta
            for portPhi = portsPhi
                doseRate3 = doseRate3 + sPass/nDoseRate(energyFlux(i)) + sourceNew/nDoseRate(sScatterN2N(1,2));
                doseRateTest3(i) = doseRateTest3(i) + sPass/nDoseRate(energyFlux(i)) + sourceNew/nDoseRate(sScatterN2N(1,2));
            end
        end
        end
        i
        
    end
    %doseRateTot = 4*(doseRate+doseRate1+doseRate3);
    
%     fid = fopen('totalDoseRate.txt', 'wt');
%     fprintf(fid, 'Dose Rate = %f\n',doseRateTot);
%     fclose(fid);
    doseRateTestTot = 4*(doseRateTest+doseRateTest1+doseRateTest3);
    doseRateMatrix = [energyFlux doseRateTestTot];
    figure
    plot(energyFlux,doseRateTestTot)
    xlswrite('doserate',doseRateMatrix)
end
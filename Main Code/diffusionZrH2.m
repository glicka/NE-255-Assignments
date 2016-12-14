function diffuseFlux = diffusionZrH2(scatterFlux,dr,mesh,yVal,thetaPrime, energy)
    flux = 0;
    for i = 1:length(dr)
        if dr(i)>0
            yInt = yVal*mesh(i,1)+tan(deg2rad(thetaPrime))*mesh(i,1);
            r_prime = neutronflux_r(yInt, thetaPrime);
            flux = flux + diff_neutronflux(scatterFlux(i),energy,r_prime);
        end
    end
    diffuseFlux = flux;
end
        
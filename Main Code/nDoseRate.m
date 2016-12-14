function Factor = nDoseRate(energy)
    [NeutronEnergy, ConversionFactor] = read_file('DoseScript.txt');
    %Neutron Energy in eV
    %Conversion Factor neutrons / (cm2*sec) per mrem/hr
    ConversionFactor = ConversionFactor * (60/(0.1*10^(-3))); %Conversion Factor is neutrons / (cm2*sec) per Gy/min
        for i = 1:length(NeutronEnergy)
            if NeutronEnergy(i) < energy && NeutronEnergy(i+1) > energy
                Conversion = ConversionFactor(i) + (ConversionFactor(i+1)-ConversionFactor(i))*(energy - NeutronEnergy(i))/(NeutronEnergy(i+1) - NeutronEnergy(i));
            elseif energy == NeutronEnergy(i)
                Conversion = ConversionFactor(i);
            end
        end   
    Factor = Conversion;
end

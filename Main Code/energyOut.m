function Energy = energyOut(data_file)
% Read the Energy values from the text file. The values extracted are energyFlux taken from serpent.
    filename = sprintf(data_file);
    fileID = fopen(filename, 'r');
    formatSpec = '%f %f %f';
    sizeA = [3 Inf];
    A = fscanf(fileID, formatSpec, sizeA);
    fclose(fileID);
    A = A';
    Energy = A(:,3);
end

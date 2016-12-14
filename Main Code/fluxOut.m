function Flux = fluxOut(data_file)
% Read the Flux values from the text file
    filename = sprintf(data_file);
    fileID = fopen(filename, 'r');
    formatSpec = '%f %f %f %f %f %f %f %f %f %f %f';
    sizeA = [12 Inf];
    A = fscanf(fileID, formatSpec, sizeA);
    fclose(fileID);
    A = A';
    Flux = A(:,11);
end

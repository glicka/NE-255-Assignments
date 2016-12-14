function [First, Second] = read_file(data_file)
% Read the data from text file.
    filename = sprintf(data_file);
    fileID = fopen(filename, 'r');
    formatSpec = '%f %f';
    sizeA = [2 Inf];
    A = fscanf(fileID, formatSpec, sizeA);
    fclose(fileID);
    A = A';
    First = A(:,1);
    Second = A(:,2);
end
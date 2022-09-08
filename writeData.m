function writeData(filename, dataVector)
    filepath = "C:\Data\Purdue\Purdue Space Program\High Altitude\Mass Estimation\PSP-HA-Structures-main\" + filename + ".txt";
    fid = fopen(filepath, 'w');
    items = size(dataVector);
    for i = 1:1:items(2)
        output = dataVector(i);
        fprintf(fid, '%f\n', output);
    end
    
    
    fclose(fid);
end
function writeData(filename, dataVector)
    filepath = "C:\Users\Amrit Arora\Desktop\Mass_Data\" + filename + ".txt";
    fid = fopen(filepath, 'w');
    items = size(dataVector);
    for i = 1:1:items(2)
        output = dataVector(i);
        fprintf(fid, '%f\n', output);
    end
    
    
    fclose(fid);
end
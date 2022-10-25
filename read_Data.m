function [output] = read_Data(filename)
    filepath = "C:\Users\Amrit Arora\Desktop\Mass_Data\" + filename + ".txt";
    fid = fopen(filepath, 'r');
    output = fscanf(fid, '%f');
end
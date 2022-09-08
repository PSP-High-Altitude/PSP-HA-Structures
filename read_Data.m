function [output] = read_Data(filename)
    filepath = "C:\Data\Purdue\Purdue Space Program\High Altitude\Mass Estimation\PSP-HA-Structures-main\" + filename + ".txt";
    fid = fopen(filepath, 'r');
    output = fscanf(fid, '%f');
end
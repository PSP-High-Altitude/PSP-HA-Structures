function writeOutputs(filename, outputs)
%WRITEOUTPUTS Summary of this function goes here
%   Detailed explanation goes here
    fid = filename + ".txt";
    writematrix(outputs, fid, 'Delimiter', 'tab');
    
end


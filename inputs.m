function inputs(MH2, MH1, D2, D1, P2, P1, acc, Drag, PL)
        mass_hist2 = cell2mat(MH2);
        mass_hist1 = cell2mat(MH1);
        MEOP2 = max(P2);
        MEOP1 = max(P1);
        max_a = max(cell2mat(acc));
        max_drag = max(cell2mat(Drag));
        
        
        writeData("test", mass_hist2);
        writeData("test2", mass_hist1);
        writeData("test3", [D2, D1, MEOP2, MEOP1, max_a, max_drag, PL]);
end


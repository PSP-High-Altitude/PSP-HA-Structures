classdef Fins < handle
    %FINS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Geometric Properties
        RC;         % Root Chord, m
        TC;         % Tip Chord, m
        SPAN;       % Span, m
        t;          % Thickness, m
        SL;         % Sweep Length, m
        FIN_COUNT;  % Number of fins, unitless
       
        % Material Properties
        rho;        % Density, kg/m^3
        
        % Design Criteria
        MASS;       % mass, kg
        CoM_y;      % CoM from fin tip, m
        CoM_x;      % CoM along fin span, m
        
    end
    
    methods
        function obj = Fins(input1,input2, input3, input4, input5, input6, input7)
            %FINS Construct an instance of this class
            %   Detailed explanation goes here
            obj.RC = input1;
            obj.TC = input2;
            obj.SPAN = input3;
            obj.t = input4;
            obj.FIN_COUNT = input5;
            obj.rho = input6;
            obj.SL = input7;
            obj.Fins_Mass(obj.RC, obj.TC, obj.SPAN, obj.t, obj.FIN_COUNT, obj.rho);
            obj.Center_of_Mass(obj.RC, obj.TC, obj.SL, obj.SPAN);
        end
        
        function Fins_Mass(obj, rootChord, tipChord, span, thickness, fincount, density)
            % finMass Calculates mass of a rocket fin
            % This function caculates the mass of a trapesodial shaped fin.
            
            area = ((rootChord + tipChord) / 2) * span;
            volume = area * thickness;
            obj.MASS = volume * density * fincount;
            
        end
        
        function Center_of_Mass(obj, RC, TC, SL, Span)
            leading_edge = @(x) -SL/Span .* x + RC;
            trailing_edge = @(x) (RC-SL-TC)/Span .* x;
            moment_fun_y = @(x, y) y;
            moment_fun_x = @(x, y) x;
            
            Moment_y = integral2(moment_fun_y, 0, Span, trailing_edge, leading_edge);
            Moment_x = integral2(moment_fun_x, 0, Span, trailing_edge, leading_edge);
            Mass = integral(leading_edge, 0, Span) - integral(trailing_edge, 0, Span);
            
            obj.CoM_y = Moment_y/Mass;
            obj.CoM_x = Moment_x/Mass;
        end
        
    end
end


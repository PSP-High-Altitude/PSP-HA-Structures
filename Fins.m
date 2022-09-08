classdef Fins < handle
    %FINS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Geometric Properties
        RC;         % Root Chord, in
        TC;         % Tip Chord, in
        SPAN;       % Span, in
        t;          % Thickness, in
        SL;         % Sweep Length, in
        FIN_COUNT;  % Number of fins, unitless
       
        % Material Properties
        rho;
        
        % Design Criteria
        MASS;
        CoM_y;
        CoM_x;
        
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
        
        %We need the fins CoM in the x-axis as well
%       function Moment_of_Inertia(obj, Mass, CoM_x, CoM_y)
        % This function uses the mass of the fins and the center of mass
        % to estimate the moment caused by the fins
        %
        % Assumptions:
        % It assumes that the fin is a point mass with a center of mass
        % (x,y)
        % Inputs:
        % Property               Variable Name           Units
        % Fin Mass                    Mass                 kG
        % Fin Center of Mass in Y     CoM_y                in
        % Fin Center of Mass in X     CoM_X                in
        %
        % Outputs(assingments):
        % Property          Variable Name           Units
        % MoI about x       moix                    lb-ft*s^2
        % MoI about y       moiy                    lb-ft*s^2
        % MoI about z       moiz                    lb-ft*s^2    
        %May not potentially need this if moment of inertia is just
        %calculated in the main function for the fins with respect to the
        %tip of the nosecone in the x, y, and z axes
      
            
        %end
        
    end
end


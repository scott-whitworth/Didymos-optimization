classdef Constants
    properties (Constant)
        MCONVERSION = 1.49597870691E11;
        SCONVERSION = 3.1557600E7;
        DEGTORADCONVERSION = pi/180
        G = 6.673848e-11*Constants.SCONVERSION^2/Constants.MCONVERSION^3;
        MS = 1.98892e30;
        ME = 5.9742e24;
        MUS = Constants.G*Constants.MS;
        MUE = Constants.G*Constants.ME;
        EARTHA = 1.49598261e11/Constants.MCONVERSION;
        EARTHE = 0.01671123; %Eccentricity of earth
        EARTHC = Constants.EARTHA*(1-Constants.EARTHE^2);
        EARTHR = 6378100 %Radius of earth (m)
        EFFBOOST =  xlsread('Input.xlsx','Operation Guidelines','b16:b16') %Percentage to boost efficiency by
        ESOI = Constants.EARTHA*(Constants.ME/Constants.MS)^.4;
        EARTHH = sqrt(Constants.EARTHC*Constants.MUS); %H value for earth (non-dimensional)

        EPSILON = 1E-6;%Small threshold
        %User Inputs
        FMIN_MAXNUMITER = 30000; %Maximum number of times fminsearch is called per time guess
        FMIN_MAXNUMEVAL = 30000;
        FMIN_FVALMAX = 1E-20; %Maximum value of F for full convergence
        FMIN_FVALMAXFAST = 1E-10; %Maximum value of F for fast convergence
        FMIN_TOLX = 1E-10; %x value tolerance for fminsearch
        FMIN_TOLFUN = 1E-10; %Function value tolerance for fminsearch
        FMIN_INTERESTINGTHRESHOLD = 1E-6; %Set interesting convergence threshold
        NUMBER_OF_INDICIES = 250;
        ODE45_RELTOL = 1e-7; %Relative tolerance of all variables for ode45
        %         ODE45_ABSTOL = [1e-6 1e-6 1e-6 1e-13]; %Absolute tolerance of each variable for ode45
        R_MIN = 0.5; %Closest the spacecraft is allowed to get to the Sun (in AU)
    end
end
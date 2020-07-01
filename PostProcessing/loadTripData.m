function [tripTime,coast_threshold,y0E,y0A,gammaCoeff,tauCoeff,coast] = loadTripData(cVector,sizeC)

%% Array offsets
        
        GAMMA_OFFSET = 1; % x[0-8] fourth order fourier for in-plane angle
        TAU_OFFSET = 8; %x[9-13] first order fourier for out-of-plane angle
        TRIPTIME_OFFSET = 14; %x[16] total duration of the trip
        COAST_OFFSET = 15; %x[17-21] second order fourier for coasting determination
        
        %% Constants
        
        tripTime=cVector(TRIPTIME_OFFSET);%Obtained from the optmization results, 15th element of the array
        coast_threshold = ones(sizeC,1)*0.5;
        
         %% Initial conditions of Earth and Asteroid
        
        % Earth
        y0E = [1.00140803662733E+00, 1.2786132931868E-01, -1.195365359889E-05,...
        -3.30528017791942E-09, 1.98791889005860E-07, -9.89458740916469E-12];
        
        % Asteroid
        y0A = [1.03524021423705E+00, 1.59070192235231E-01, -5.54192740243213E-02,...
        -2.53512430430384E-08, 2.27305994342265E-07, 7.29854417815369E-09];

        
        %% fourier components
        
        gammaCoeff= cVector(GAMMA_OFFSET:GAMMA_OFFSET+6);
        tauCoeff=cVector(TAU_OFFSET:TAU_OFFSET+4);
        coast=cVector(COAST_OFFSET:COAST_OFFSET+4);
end


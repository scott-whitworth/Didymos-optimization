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
        y0E = [1.00021392223428, 0.199470650149394, -1.54878511585620e-05,...
            -3.32034068725821e-09, 1.99029138292504e-07, -9.71518257891386e-12];
        
        % Asteroid
        y0A = [1.02696822710421, 0.238839574416454, -0.0526614832914496,...
            -2.05295246185041e-08, 2.29132593453064e-07, 8.00663905822009e-09];

        
        %% fourier components
        
        gammaCoeff= cVector(GAMMA_OFFSET:GAMMA_OFFSET+6);
        tauCoeff=cVector(TAU_OFFSET:TAU_OFFSET+4);
        coast=cVector(COAST_OFFSET:COAST_OFFSET+4);
end


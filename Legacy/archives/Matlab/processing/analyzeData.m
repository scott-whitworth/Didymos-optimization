function [] = analyzeData(runLength)
format longE
%% Load file
runDataExport=zeros(3,runLength);
cSize = 20;
cVectorExport=zeros(cSize-1,runLength);
for i=0:runLength
    filenameO=join(['final-optimization',num2str(i-1),'.bin']);
    filenameT=join(['orbitalMotion-accel',num2str(i-1),'.bin']);
    if isfile(filenameO)
        
        fileOptimized = fopen(filenameO);
        cVector = fread(fileOptimized,cSize,'double');
        fileY = fopen(filenameT);
        sizeC=cVector(end);
        cR = fread(fileY,[11 sizeC],'double');
        
        %% Load mission data
        
        [tripTime,coast_threshold,y0E,y0A,gammaCoeff,tauCoeff,coast] = loadTripData(cVector,sizeC);
                
        %% comparison of final velocity
        
        AU = 1.49597870691e11; % m/AU
        errorPos = sqrt((y0A(1,1)-cR(1,end))^2+(y0A(1,2)-cR(2,end))^2+(y0A(1,3)-cR(3,end))^2)*AU;
        finalVel = sqrt((y0A(1,4)-cR(4,end))^2+(y0A(1,5)-cR(5,end))^2+(y0A(1,6)-cR(6,end))^2)*AU;
        %missionFinVel = 4.399378072e-08*AU;
        %display(finalVel);
        fuel=cR(11,end-1);
        %display(fuel)
        
        %plotData(cR,y0A,y0E,sizeC,timeFinal,tripTime,coast,coast_threshold)
        
        %% clearing variables
        
        fclose('all');
        
        
        %% make folder for new data
        
        dataName=join([num2str(finalVel),'-',num2str(fuel),'-',num2str(tripTime)]);
        newFolder=join([pwd,'\results','\',dataName]);

        if ~exist(newFolder, 'dir')
            mkdir(newFolder);
        end
        
        % Move the file to its folder
        movefile(filenameO,newFolder)
        movefile(filenameT,newFolder)
        
        % Record solution stats
        runDataExport(1,i)=finalVel;
        runDataExport(2,i)=fuel;
        runDataExport(3,i)=tripTime;
        
        cVectorExport(:,i)=cVector(1:end-1);
        
        
    else
        %disp('No file found')
    end
    
    clear a accelR accelTheta accelX accelY accelZ ALPHA_OFFSET BETA_OFFSET GAMMA_OFFSET TAU_OFFSET TRIPTIME_OFFSET ...
            aX aY aZ b  cVX cVY cVZ cX cY cZ eX eY eZ fileC fileY gamma gammaCoeff options radStep tA tau tauCoeff tE ...
            timeFinal tspan y0A y0E AU
end

clear cSize runLength

% clear runDataExport
runDataExport( :, all(~runDataExport,1) ) = [];
runDataExport=runDataExport.';
cVectorExport( :, all(~cVectorExport,1) ) = [];
cVectorExport=cVectorExport.';
% write data to export to the application needed
dlmwrite('csvSummary.txt',runDataExport,'delimiter',',','precision','%.14d')
%dlmwrite('cGuess.txt',cVectorExport,'delimiter',',','precision','%.14d')
fileVectorOptimized = fopen('optimizedVector.bin','w');
fwrite(fileVectorOptimized,cVectorExport,'double');
fclose(fileVectorOptimized);
end


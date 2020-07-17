function [] = filePlot(rank,seed)

filenameO=join(['bin/finalOptimization-',num2str(rank),'-',num2str(seed),'.bin']);
filenameT=join(['bin/orbitalMotion-',num2str(rank),'-',num2str(seed),'.bin']);
fileC = fopen(filenameO);
cVector = fread(fileC,Inf,'double');
fileY = fopen(filenameT);
sizeC=cVector(end)+1;
cR = fread(fileY,[11, sizeC],'double');
[tripTime,coast_threshold,y0E,y0A,gammaCoeff,tauCoeff,coast] = loadTripData(cVector);
plotData(cR,y0A,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff)
fclose('all');
end


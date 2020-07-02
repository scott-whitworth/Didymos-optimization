function [] = filePlot(seed, rank)
cSize=20;

filenameC=join(['configuration-',num2str(seed),'.bin']);
filenameO=join(['finalOptimization-',num2str(seed),'-',num2str(rank),'.bin']);
filenameT=join(['orbitalMotion-',num2str(seed),'-',num2str(rank),'.bin']);
fileConfig = fopen(filenameC);
config = fread(fileConfig,Inf,'double');
fileC = fopen(filenameO);
cVector = fread(fileC,cSize,'double');
fileY = fopen(filenameT);
sizeC=cVector(end)+1;
cR = fread(fileY,[11, sizeC],'double');
[tripTime,coast_threshold,y0E,y0A,gammaCoeff,tauCoeff,coast] = loadTripData(cVector,sizeC,config);
plotData(cR,y0A,y0E,sizeC,tripTime,coast,coast_threshold,gammaCoeff,tauCoeff)
fclose('all');
end


function [] = filePlotCompare(rank1,seed1,rank2,seed2)

fileC = fopen(join(['bin/finalOptimization-',num2str(rank1),'-',num2str(seed1),'.bin']));
cVector = fread(fileC,Inf,'double');
fileY = fopen(join(['bin/orbitalMotion-',num2str(rank1),'-',num2str(seed1),'.bin']));
sizeC=cVector(end)+1;
fileD = fopen(join(['bin/finalOptimization-',num2str(rank2),'-',num2str(seed2),'.bin']));
dVector = fread(fileD,Inf,'double');
fileZ = fopen(join(['bin/orbitalMotion-',num2str(rank2),'-',num2str(seed2),'.bin']));
sizeD = dVector(end)+1;
cR = fread(fileY,[11, sizeC],'double');
dR = fread(fileZ,[11, sizeD],'double');
[tripTime1,coast_threshold1,y0E,y0A,gammaCoeff1,tauCoeff1,coast1] = loadTripData(cVector);
[tripTime2,coast_threshold2,y0E,y0A,gammaCoeff2,tauCoeff2,coast2] = loadTripData(dVector);
plotDataCompare(cR,y0A,y0E,sizeC,tripTime1,coast1,coast_threshold1,gammaCoeff1,tauCoeff1,dR,sizeD,tripTime2,coast2,coast_threshold2,gammaCoeff2,tauCoeff2)
fclose('all');
end

function [] = filePlotCompare(seed1,seed2)
size=36;

fileC = fopen(join(['finalOptimization-',num2str(seed1),'.bin']));
cVector = fread(fileC,size,'double');
fileY = fopen(join(['orbitalMotion-',num2str(seed1),'.bin']));
sizeC=cVector(end)+1;
fileD = fopen(join(['finalOptimization-',num2str(seed2),'.bin']));
dVector = fread(fileD,size,'double');
fileZ = fopen(join(['orbitalMotion-',num2str(seed2),'.bin']));
sizeD = dVector(end)+1;
cR = fread(fileY,[11, sizeC],'double');
dR = fread(fileZ,[11, sizeD],'double');
[tripTime1,coast_threshold1,y0E,y0A,gammaCoeff1,tauCoeff1,coast1] = loadTripData(cVector,sizeC);
[tripTime2,coast_threshold2,y0E,y0A,gammaCoeff2,tauCoeff2,coast2] = loadTripData(dVector,sizeD);
plotDataCompare(cR,y0A,y0E,sizeC,tripTime1,coast1,coast_threshold1,gammaCoeff1,tauCoeff1,dR,sizeD,tripTime2,coast2,coast_threshold2,gammaCoeff2,tauCoeff2)
fclose('all');
end

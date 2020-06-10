clear

%% To analyze several files
%analyzeData(300);

%% To plot data

dataNum=1;
filePlot(dataNum)

numberSolutions=14;
[file,error]=fopen('optimizedVector.bin');
a = fread(file,[numberSolutions 19],'double');
fclose(file);

hold on
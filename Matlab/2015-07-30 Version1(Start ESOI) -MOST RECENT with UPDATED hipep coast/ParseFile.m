%PARSEFILE Reads in file from JPL Horizons and returns vector of x and y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function reads in a file containing x and y coordinates for the
%asteroid and Earth and returns vectors of those coordinates.  It is
%designed to read a file of a specific format. It is important to ensure
%that the tables given as inputs cover the same time span with the same
%interval. This program requires the position vectors in units of km and
%velocity vectors in km/sec. Specific file names are described in the
%comments below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by GetCloseApproach.
function [Xval,Yval,AsteroidEccentricity,AsteroidSemimajor,AsteroidPeriod,zAxisFraction] = ParseFile(body,astDes)
%
%Tells the program which file to read from
if strcmpi(body,'EarthPosition')
    %If reading Earth position values, the file name must be
    %'EarthPosition.txt'
    fid = fopen(strcat('EarthData/',astDes,'Pos.txt'));
elseif strcmpi(body,'EarthVelocity')
    %If reading Earth velocity values, the file name must be
    %'EarthVelocity.txt'
    fid = fopen(strcat('EarthData/',astDes,'Vel.txt'));
elseif strcmpi(body,'AsteroidPosition')
    %If reading Asteroid position values, the file name must be
    %'AsteroidPosition.txt'
    fid = fopen(strcat('AsteroidData/',astDes,'Pos.txt'));
elseif strcmpi(body,'AsteroidVelocity')
    %If reading Asteroid velocity values, the file name must be
    %'AsteroidVelocity.txt'
    fid = fopen(strcat('AsteroidData/',astDes,'Vel.txt'));
end

%Reads the text file one line at a time, keeping the whole line as a string
%'temp' is a dummy cell array--it is not actually usable
temp = textscan(fid, '%s','delimiter','\n');
fclose(fid);

%Reads values from 'temp' into a usable cell array 'RawInput' containing
%each line of the original data file as a string
RawInput = temp{1};

%Eccentricity and semimajor axis are both non-dimensional in the file
AsteroidEccentricity = 0;
AsteroidSemimajor = 0;
AsteroidPeriod = 0;
%Find eccentricity and semimajor axis values
if strcmpi(body,'AsteroidPosition')
    for i = 1:1:size(RawInput,1)
        if(findstr(RawInput{i},'PER=') > 0)
            [~,remainingTerm] = strtok(RawInput{i});
            [secondTerm] = strtok(remainingTerm);
            [thirdTerm] = strtok(secondTerm);
            [fourthTerm] = strtok(thirdTerm);
            AsteroidPeriod = str2double(fourthTerm);
            break;
        end
    end
    
    for i = 1:1:size(RawInput,1)
        if(findstr(RawInput{i},'EC=') > 0)
            %i is the line at which eccentricity is found in txt file
            %strtok removes first term from string and puts it into
            %firstterm, keep the rest of the string in remainingTerm
            [~,remainingTerm] = strtok(RawInput{i});
            [secondTerm] = strtok(remainingTerm);
            [thirdTerm] = strtok(secondTerm);
            AsteroidEccentricity = str2double(thirdTerm);
            break;
        end
    end
    
    for i = 1:1:size(RawInput,1)
        if(findstr(RawInput{i},'A=') > 0)
            [~,remainingTerm] = strtok(RawInput{i});
            [secondTerm] = strtok(remainingTerm);
            AsteroidSemimajor = str2double(secondTerm);
            break;
        end
    end
    for i = 1:1:size(RawInput,1)
        if(findstr(RawInput{i},'ZAX=') > 0)
            [~,remainingTerm] = strtok(RawInput{i});
            [secondTerm] = strtok(remainingTerm);
            [thirdTerm] = strtok(secondTerm);
            [fourthTerm] = strtok(thirdTerm);
            zAxisFraction = str2double(fourthTerm);
            break;
        end
    end
end

%Runs through the unformatted data file and identifies the index at which
%the data starts. 'HeadEnd' is a logical array containing a '1' value at
%the index where the '$$SOE' line is located.
%The end of the header is indicated by the line '$$SOE' in the data file.
HeadEnd = strcmp('$$SOE',RawInput);

%Loop through arrays to find where each data begins
for k = 1:1:size(HeadEnd,1);
    if HeadEnd(k) == 1
        %The data begins the line AFTER the '$$SOE' character, but the
        %first line of data is text, not useful numbers.
        DataStart = k+2;
    end
end

%Runs through the unformatted data file to find the index at which the data
%ends, which occurs right before a line containing '$$EOE'. This works much
%the same as the section above.
FootStart = strcmp('$$EOE',RawInput);

for n = 1:1:size(FootStart,1);
    if FootStart(n) == 1
        %The data ends the line BEFORE the '$$EOE' character
        DataEnd = n-1;
    end
end

%Creates a cell of the data values only
%NOTE: The expression (DataEnd-DataStart)/2+1 gives the number of lines of
%usable data contained in the data file. There is a date which is discarded
%every other line, and the data set begins and ends with usable data 
%(hence the +1)
data = cell((DataEnd-DataStart)/2+1,1);

%This is an index indicator for use in the next loop
%I wish this were nicer, but I don't know what else I can do
m = 1;

%Reads coordinates from cell array
for i = DataStart:2:DataEnd
    data{m} = str2num(RawInput{i});
    m = m+1;
end

%Reads the current XYZ coordinates from the cell array
for j = 1:1:size(data,1)
    
    %Reads the line of data into a vector containing X,Y,Z data in double format
    XYZvector = data{j};
    
    %Creates vectors containing all x and y values. These will be either
    %positions or velocities, depending on which file is being read from.
    %X value is the first value in the line from the data file
    Xval(j) = XYZvector(1);
    %Y value is the second value
    Yval(j) = XYZvector(2);
    %NOTE: Z values are in the third index. Since this program works only
    %in two dimensions, the Z components are simply discarded.
end

end
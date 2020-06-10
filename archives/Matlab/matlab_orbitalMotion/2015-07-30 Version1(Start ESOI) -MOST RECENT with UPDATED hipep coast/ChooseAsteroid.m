function [astDes,astChoice] = ChooseAsteroid(dataEntry)
%CHOOSEASTEROID Allows for user input of asteroid.

%Open list of available thruster data and create cell array inputAsteroids
%to hold the data.
savedAsteroids = fopen('SavedAsteroids.txt');
temp = textscan(savedAsteroids,'%s','delimiter','/r/n');
fclose(savedAsteroids);
inputAsteroids = temp{1};
fileSize = size(inputAsteroids,1);

%Output possible asteroid choices
fprintf('\n');
if dataEntry == 0
fprintf('Input asteroid ID to be evaluated:\n');
for i = 1:1:fileSize
    fprintf('\t');
    fprintf('%i',i);
    fprintf(' for ');
    fprintf('%s\n',inputAsteroids{i});
    if i == fileSize
        fprintf('\n');
        fprintf('\t');
        fprintf('%i',i+1);
        fprintf(' to ');
        fprintf('enter new asteroid\n');
        
        fprintf('\t');
        fprintf('%i',i+2);
        fprintf(' to ');
        fprintf('delete asteroid\n');
    end
end
end
%If there are no asteroids in the data file, prompt the user to input a new
%asteroid.
if fileSize == 0
    fprintf('Input asteroid ID to be evaluated:\n');
    fprintf('\t');
    fprintf('1');
    fprintf(' to ');
    fprintf('enter new asteroid\n');
end

astChoice = fileSize + 3;
newFlag = 0;

%If the user chooses to enter a new asteroid or delete an asteroid,
%continue to output choices.
while astChoice + newFlag > fileSize
    newFlag = 0;
    if dataEntry == 0
    astChoice = input(' : ');
    else
       astChoice = xlsread('Input.xlsx','Base Inputs','b5:b5');
    end
    %If the user chooses an asteroid, set astDes to that string.
    if astChoice <= fileSize
        astDes = inputAsteroids{astChoice};
        %If the user chooses to input a new asteroid, open the data file and
        %write the new asteroid name to the file.
    elseif astChoice == fileSize + 1
        newAsteroid = input('Input asteroid designation (case and space sensitive, i.e. 1996XB27): ','s');
        saveNewAsteroid = fopen('SavedAsteroids.txt','a');
        fprintf(saveNewAsteroid,'%s\r\n',newAsteroid);
        fclose(saveNewAsteroid);
        newFlag = 1;
        fprintf('\nPlease run HORIZONS Data Retrieval Tool v2 before running running the new asteroid.\n\n');
        %If the user chooses to delete an asteroid, prompt for an asteroid to
        %delete, delete that entry from the data file, and rewrite the
        %remaining asteroid choices to the file.
    elseif astChoice == fileSize + 2
        deleteAsteroid = input('Input asteroid ID to delete: ');
        inputAsteroids{deleteAsteroid} = [];
        saveDeletedAsteroid = fopen('SavedAsteroids.txt','w');
        for j = 1:fileSize
            fprintf(saveDeletedAsteroid,'%s\r\n',inputAsteroids{j});
        end
        fclose(saveDeletedAsteroid);
    end
    
    %If the user chose to delete an asteroid, prompt the user for another
    %choice.
    if astChoice > fileSize
        
        %Update the array inputAsteroids containing the available asteroid
        %choices.
        savedAsteroids = fopen('SavedAsteroids.txt');
        temp = textscan(savedAsteroids,'%s','delimiter','/r/n');
        fclose(savedAsteroids);
        inputAsteroids = temp{1};
        fileSize = size(inputAsteroids,1);
        
        %Output asteroid choices.
        fprintf('\n');
        fprintf('Input asteroid ID to be evaluated:\n');
        for i = 1:1:fileSize
            fprintf('\t');
            fprintf('%i',i);
            fprintf(' for ');
            fprintf('%s\n',inputAsteroids{i});
            if i == fileSize
                fprintf('\n');
                fprintf('\t');
                fprintf('%i',i+1);
                fprintf(' to ');
                fprintf('enter new asteroid\n');
                
                fprintf('\t');
                fprintf('%i',i+2);
                fprintf(' to ');
                fprintf('delete asteroid\n');
            end
        end
    end
end


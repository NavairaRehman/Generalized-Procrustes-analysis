% Input and Output directories
inputFolder = '/Users/navairarehman/Desktop/rsme2';
outputFolder = '/Users/navairarehman/Documents/FYP/3d preprocess/Procrustes/matlab_test/';

% Get list of .obj files in the input folder
files = dir(fullfile(inputFolder, '*.obj'));

% Initialize cell array for storing vertex sets for GPA
allVertices = {};

% Step 1: Read all OBJ files and extract vertices and faces
for i = 1:length(files)
    % Read the OBJ file
    filePath = fullfile(inputFolder, files(i).name);
    [vertices, faces] = readObj(filePath);
    
    % Store vertices for GPA
    allVertices{i} = vertices;
end

% Step 2: Perform Generalized Procrustes Analysis (GPA)
[transformedVertices, ~, ~] = generalizedProcrustes(allVertices);

% Step 3: Save the transformed OBJ files
for i = 1:length(files)
    % Get file name
    fileName = files(i).name;
    
    % Get transformed vertices
    vertices = transformedVertices{i};
    
    % Get original faces
    filePath = fullfile(inputFolder, files(i).name);
    [~, faces] = readObj(filePath);
    
    % Save the transformed OBJ file
    saveObj(vertices, faces, fullfile(outputFolder, fileName));
end

disp('Processing complete. Transformed files saved.');


function [vertices, faces] = readObj(filePath)
    vertices = [];
    faces = {};
    
    fid = fopen(filePath, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, 'v ')
            % Extract vertex coordinates
            vertices = [vertices; sscanf(line(3:end), '%f %f %f')'];
        elseif startsWith(line, 'f ')
            % Retain face line as is
            faces{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);
end

function [transformedVertices, meanShape, scale] = generalizedProcrustes(allVertices)
    nFiles = length(allVertices);
    nPoints = size(allVertices{1}, 1);
    
    % Initialize mean shape
    meanShape = allVertices{1};
    
    % Iteratively update mean shape
    for iter = 1:10
        alignedShapes = zeros(nFiles, nPoints, 3);
        for i = 1:nFiles
            % Align to current mean shape
            [d, Z, transform] = procrustes(meanShape, allVertices{i});
            alignedShapes(i, :, :) = Z;
        end
        % Update mean shape
        meanShape = squeeze(mean(alignedShapes, 1));
    end
    
    % Final alignment
    transformedVertices = cell(nFiles, 1);
    for i = 1:nFiles
        [d, Z, transform] = procrustes(meanShape, allVertices{i});
        transformedVertices{i} = Z;
    end
    
    scale = sqrt(sum(meanShape(:).^2)); % Scale factor for normalization
end


function saveObj(vertices, faces, outputPath)
    fid = fopen(outputPath, 'w');
    
    % Write vertices
    for i = 1:size(vertices, 1)
        fprintf(fid, 'v %.6f %.6f %.6f\n', vertices(i, 1), vertices(i, 2), vertices(i, 3));
    end
    
    % Write faces
    for i = 1:length(faces)
        fprintf(fid, '%s\n', faces{i});
    end
    
    fclose(fid);
end

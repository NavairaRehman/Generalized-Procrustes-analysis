
%%  This piece of code employs GPA for face alignment : developed by Pulak Purkait 
%   The details of the algorithms can be found in 
%   https://graphics.stanford.edu/courses/cs164-09-spring/Handouts/paper_shape_spaces_imm403.pdf
%   Email - pulak.isi@gmail.com 

%% Read datasets  
data_path = '/Users/navairarehman/Documents/FYP/3d preprocess/Procrustes/300w_cropped/01_Indoor/'; 
filelist = dir(data_path);

%Calculating the # of datasample (images and their corresponding pts files)
datasize = numel(filelist)/2 -1; 

% Note that the images are of different size 
% the following lines if image has the same dimensions. 
% sz = size(imread([data_path, filelist(3).name])); 
% images = zeros([datasize, sz]); 

flag_display = 0; % You can choose this to be 0 if you dont want to display 
images = cell(1, datasize); %using a cell array for images to store different dimensions of image data. 

%reading a pts file to get the no. of points/landmarks. They are stored in
%the header of the pts files of the test dataset they have used.
%68 in this case
fileID = fopen([data_path, filelist(4).name], 'r'); 
fscanf(fileID, '%s', 3); 
no_points = fscanf(fileID, '%d\n', 1); 
fclose(fileID);

%3 dim array to store x and y coordinates of all the points of each
%datapoint in the dataset.
points = zeros([datasize, 2, no_points]); 

%readign the image and pts file and populating the points array and image
%cell array
for i = 1:datasize 

    images{i} = imread([data_path, filelist(2*i+1).name]); 

    fileID = fopen([data_path, filelist(2*i+2).name], 'r'); 
    fscanf(fileID, '%s', 3); 
    no_points = fscanf(fileID, '%d\n', 1); 

    fscanf(fileID, '%s', 1); 
    points(i, :, :) = reshape(fscanf(fileID, '%f', 2*no_points), [2, no_points]); 
    disp(i/datasize); 
    disp("done")

% Plot : To check if the readings are correct 
    if flag_display
        imshow( images{i} ); 
        hold on; 
        plot(squeeze(points(i, 1, :)),  squeeze(points(i, 2, :)), 'r*'); 
        hold off; 
        %pause; %comment this out. It enables viewing the landmarks on each photo if flag_display in on. 
    end

end 

% save('points.mat', 'points'); 
% Procrustes analysis https://en.wikipedia.org/wiki/Procrustes_analysis 
% load('points.mat'); 
% Translation (so that the mean of the points is cetered at the origin)

points_mean = repmat(mean(points, 3), [1, 1, no_points]); 
points = points - points_mean; 

%% Scaling (to have the same size by adjusting for differences distance)

points_scale = sqrt(sum(sum(points.^2, 3), 2)/no_points); 
points = points./repmat(points_scale, [1, 2, no_points]); 

%% Rotation (to minimize the distance and difference between the points)

theta = atan(squeeze(sum(points(:, 1, :).*repmat(points(1, 2, :), [datasize, 1, 1]) - points(1, 2, :).*repmat(points(1, 1, :), [datasize, 1, 1]), 3))...
    ./squeeze(sum(points(:, 1, :).*repmat(points(1, 1, :), [datasize, 1, 1]) + points(1, 2, :).*repmat(points(1, 2, :), [datasize, 1, 1]), 3))); 

R1 = [cos(theta), -sin(theta)]; 
R2 = [sin(theta), cos(theta)]; 

points_transformed = cat(2, sum(points.*repmat(R1, [1, 1, no_points]), 2), sum(points.*repmat(R2, [1, 1, no_points]), 2)); 
disp("procrustes done")
disp(points_transformed)

% Checking again if the things are going well 
%displaying the transformed points after translation, scaling and rotation.
if flag_display
    for i = 1:datasize 
        plot(squeeze(points_transformed(i, 1, :)), - squeeze(points_transformed(i, 2, :)), 'r*'); 
        hold off; 
        %pause;  
    end
end

%% Generalized Procrustes analysis 
% Compute optimal mean face C_hat
% Initialize the initial centroid C as the first facial image 

% load('points.mat'); 
% Translation 
disp("Performing GPA")
no_points = size(points, 3);   
datasize = size(points, 1);   

points_mean = mean(points, 3); 
points_mean = repmat(points_mean, [1, 1, no_points]); 

points = points - points_mean; 


N = 1000;                                % Maximum number of iterations 
epsln = 0.01;                            % Quality Parameter 
C_hat = squeeze(points(1, :, :));        %initial estimate of the mean shape. Taken as the first image here.
Error_plot = zeros(N, 1); 
Error_plot(1) = 1; 
iter = 1; 
disp("going into GPA loop")
while (Error_plot(iter) > epsln && iter < N) 
    C = C_hat; 
    % Compute similarity transformations between the images and mean shape 
    C_hat = zeros(2, no_points); 
    for i = 1:datasize
        X = squeeze(points(i, :, :)); 
        s = sqrt(sum(C(:).^2) / sum(X(:).^2)); % This will give us the relative scale 
        M = C*X';
        R = M*(M'*M)^(-0.5);                   % 
        % Fitting a similarity transform between the target shape and
        % original shape 
        C_hat = C_hat + s*R*X; 
    end
    C_hat = C_hat/datasize; 
    C_hat = C_hat - repmat(mean(C_hat, 2), [1, no_points]); 
    iter = iter + 1; 
    Error_plot(iter) = norm(C - C_hat); 
    disp(Error_plot(iter)); 
end
Error_plot(iter+1:end) = []; 
plot(Error_plot); % Ploting Error values 

%% Translate the points with optimal centre 

points_transformed = points; 

for i = 1:datasize
    X = squeeze(points(i, :, :)); 
    s = sqrt(sum(C_hat(:).^2) / sum(X(:).^2)); %scaling factor caluclated for each image  
    M = C_hat*X';
    R = M*(M'*M)^(-0.5); %rotation matrix for each image
    points_transformed(i, :, :) = s*R*squeeze(points(i, :, :)); 
end

%% Checking again if the things are going well 
if flag_display
    for i = 1:datasize 
        plot(squeeze(points_transformed(i, 1, :)), - squeeze(points_transformed(i, 2, :)), 'r*'); hold on; 
        plot(C_hat(1, :), - C_hat(2, :), 'g*'); 
        axis([-0.15 0.15 -0.15 0.15]); 
        hold off; 
        %pause; 
    end
end
% 
% %% Principle component Analysis of shapes 
% X = reshape(points_transformed, [datasize, 2*no_points]); 
% X = bsxfun(@minus,X,mean(X));
% % [coeff,score,latent] = pca(X); 
% [V, D] = eig(X'*X); 
% D=diag(D);
% [D_sorted, id]=sort(D, 'descend'); % sort according to the magnitude of eigenvalues 
% D_sum = cumsum(D_sorted)/sum(D_sorted); 
% id_c = find(D_sum < 0.98);        % keep 98% data 
% 
% M = V(:, id(id_c));               % Eigenvectors corresponding to large eigen values 
% Y = X*M;                          % Data projected on the eigenspaces 
% 
% eigenshapes = reshape(M, [numel(id_c), 2, no_points]); % Computing Eigen Shapes 
% 
% % Plotting the Eigen Shapes 
% if flag_display
%     for i = 1:datasize 
%         plot(squeeze(eigenshapes(i, 1, :)), - squeeze(eigenshapes(i, 2, :)), 'r*'); hold on; 
%         plot(C_hat(1, :), - C_hat(2, :), 'g*'); 
%         axis([-0.15 0.15 -0.15 0.15]); 
%         hold off; 
%         pause; 
%     end
% end
% X_bar = Y*M';                    % Back to the original shape 
% 
% disp(norm(X - X_bar)); 

% saving the points in a pts file 
% Define output path where the transformed pts files will be saved
output_path = '/Users/navairarehman/Documents/FYP/3d preprocess/Procrustes/300w_cropped/aligned_pts/';

% Loop through all the images and save the transformed points into new .pts files
for i = 1:datasize
    % Create the new .pts file name (you can modify this as per your requirement)
    [~, filename, ~] = fileparts(filelist(2*i+2).name); % Extract filename without extension
    output_filename = [filename, '_aligned.pts']; % Create new filename for aligned points
    
    % Open the new .pts file for writing
    fileID = fopen([output_path, output_filename], 'w');
    
    % Write the header
    fprintf(fileID, 'version: 1\n');
    fprintf(fileID, 'n_points: %d\n', no_points);
    fprintf(fileID, 'point\n');
    
    % Write the transformed points to the file
    transformed_points = squeeze(points_transformed(i, :, :)); % Extract the transformed points for the current image
    fprintf(fileID, '%f %f\n', transformed_points); % Write the points as x y coordinates
    
    % Close the file
    fclose(fileID);
    
    disp(['Saved aligned points for image ', num2str(i)]);
end

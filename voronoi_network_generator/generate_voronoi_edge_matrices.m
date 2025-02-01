close all
clear all
clc

format long;
% domain size in x and y directions (I think of these as if in microns)
bound_x = 1105.0;  
bound_y = 1105.0 + (2.0*85.0);
% bound_x = 2000;  
% bound_y = 2000;
% polybdn_voronoi() function below creates voronoi diagram with edges
% restricted to our domain only (normal Voronoi diagrams sticks out of a given domain quite a bit)
BdM=[0,0;bound_x,0;bound_x,bound_y;0,bound_y;0,0];

% Define the minimum edge length
min_length = 7.39;

% Set the number of seed points for the Voronoi diagram
SeedPointsN = 60;
start_Sample = 959;
SampleN = 1250; % number of samples to generate

% Set the random seed
base_seed = 48;
% redo_selections = [134, 227, 335, 355, 369, 376, 392, 395, 408, 419, 429, 446, 468, 536, 546, 572, 587, 589, 591, 621, 634, 639, 641, 648, 654, 658, 688, 713, 714, 747, 754, 766, 773, 776, 790, 807, 821, 822, 843, 846, 856, 900, 909, 913, 942, 950, 953, 958, 959, 961, 963, 973, 977, 986, 997];
% for k = redo_selections
for k=start_Sample:SampleN
    
    % Define the number of attempts
    max_attempts = 2000; % Adjust this number based on how many attempts we want
    attempt = 0;
    
    % Initialize variables for the while loop
    valid_edges_found = false;
    
    while ~valid_edges_found && attempt < max_attempts
        
        try
            attempt = attempt + 1;

            % Combine the base seed with the loop indices
            seed = base_seed + 10000 * (k - 1) + attempt;
            
            % Use modulo to keep the seed within a reasonable range
            % This example uses a prime number (2^31 - 1) as the modulus
            current_seed = mod(seed, 2147483647);
            rng(current_seed);
            % Generate the first random number
            seed_x = randi(1000000);  % Adjust the range as needed
           
            % The second random number is always 1 more than the first
            seed_y = seed_x + 1;
            % current_seed_factor = random_number + attempt;
            % current_seed = 2*(current_seed_factor-1);
            % rng(current_seed); % Use a different seed for each attempt
            
            disp('Attempt:');
            disp(attempt);
            disp(k);
            disp(seed_x);
            disp(seed_y);
            % disp(valid_edges_found);
    
            % Pick SeedPointsN (number) of random points in the square domain
            % rng(seed_x);  % Use the seed to control random number generation
            x = bound_x * rand(1, SeedPointsN);  % Generate random x-coordinates
            % rng(seed_y);  % Use another seed to control the y-coordinates separately
            y = bound_y * rand(1, SeedPointsN);  % Generate random y-coordinates
            % x = bound_x*gallery('uniformdata',[1 SeedPointsN],seed_x);
            % y = bound_y*gallery('uniformdata',[1 SeedPointsN],seed_y);
            P=[x',y'];

            % This generates a Voronoi diagram constrained to a given domain
            [vornb,vorvx,Aaug,baug] = polybnd_voronoi(P,BdM);
            A = [];
            
            for i = 1:size(P,1)
                 for j=1:(length(vorvx{i}(:,1))-1)
                     if (round(vorvx{i}(j,1),10)==0 & round(vorvx{i}(j+1,1),10)==0)|(round(vorvx{i}(j,1),10)== bound_x & round(vorvx{i}(j+1,1),10)==bound_x)|(round(vorvx{i}(j,2),10)==0 & round(vorvx{i}(j+1,2),10)==0)|(round(vorvx{i}(j,2),10)==bound_y & round(vorvx{i}(j+1,2),10)==bound_y)
                     else 
                         A=round([A;vorvx{i}(j,1),vorvx{i}(j,2),vorvx{i}(j+1,1),vorvx{i}(j+1,2)],10);
                     end
                 end
            end
        
            cond = (A(:,1)>A(:,3))|(A(:,1)==A(:,3) & A(:,2)>A(:,4));
            A(cond,:) = [];
            
            [~,idx] = sort(A(:,1));
            Adj = A(idx,:);
            
            AdjWithoutSides=Adj;
            cond2 = (round(AdjWithoutSides(:,2),10)==0.0)|(round(AdjWithoutSides(:,2))==bound_y)|(round(AdjWithoutSides(:,4),10)==0.0)|(round(AdjWithoutSides(:,4))==bound_y);
            AdjWithoutSides(cond2,:) = [];
            
            % Condition 1: Find edges where x1 == 0
            condition1 = AdjWithoutSides(:, 1) == 0;
            if any(condition1)
                % Extract edges where x1 == 0
                subset = AdjWithoutSides(condition1, :);
                
                % Compute the average y1 of these edges
                avg_y1 = mean(subset(:, 2));
                
                % Find the edge where y1 is closest to the average y1
                [~, idx_closest] = min(abs(subset(:, 2) - avg_y1));
                edge_to_keep = subset(idx_closest, :);
                
                % Remove all edges where x1 == 0 and keep only the selected edge
                AdjWithoutSides = AdjWithoutSides(~condition1, :);
                AdjWithoutSides = [edge_to_keep; AdjWithoutSides];
            end

            % Condition 2: Find edges where x2 ~= bound_x
            condition2 = AdjWithoutSides(:, 3) == bound_x;
            if any(condition2)
                % Extract edges where x2 ~= bound_x
                subset = AdjWithoutSides(condition2, :);
                
                % Compute the average y2 of these edges
                avg_y2 = mean(subset(:, 4));
                
                % Find the edge where y2 is closest to the average y2
                [~, idx_closest] = min(abs(subset(:, 4) - avg_y2));
                edge_to_keep = subset(idx_closest, :);
                
                % Remove all edges where x2 ~= bound_x and keep only the selected edge
                AdjWithoutSides = AdjWithoutSides(~condition2, :);
                AdjWithoutSides = [AdjWithoutSides; edge_to_keep]; % Ensure edge_to_keep is at the bottom
            end

            % Calculate the lengths of the remaining edges
            filtered_edges = AdjWithoutSides;
            x1 = filtered_edges(:, 1);
            y1 = filtered_edges(:, 2);
            x2 = filtered_edges(:, 3);
            y2 = filtered_edges(:, 4);
            edge_lengths = sqrt((x2 - x1).^2 + (y2 - y1).^2);
            disp(min(edge_lengths));
            % Check if ALL edges meet the minimum length requirement
            if all(edge_lengths >= min_length)
                valid_edges = filtered_edges;
                valid_edges_found = true; % Exit loop if all edges are valid
            end
        catch ME
            fprintf('Error in processing at Sample %d: %s\n', k, ME.message);
            continue; 
        % If any edge is smaller than min_length, the loop will continue
        % and generate a new distribution
        end
    end




    if valid_edges_found
        rootdir = '/home/narain/Desktop/Scripts/voronoi_network_generator/output';
        relativefolder = append(num2str(SeedPointsN), 'SeedPoints');
        fullfolderpath = fullfile(rootdir, relativefolder);
        
        % Create the directory if it does not exist
        if ~exist(fullfolderpath, 'dir')
            mkdir(fullfolderpath);
        end
        
        % Open the file for writing
        EdgesMatrixFile = fopen(fullfile(fullfolderpath, append('EdgesMatrixSampleNumber', num2str(k), '.txt')), 'wt');
        
        for ii = 1:size(valid_edges,1)
            fprintf(EdgesMatrixFile,'%g\t',valid_edges(ii,:));
            fprintf(EdgesMatrixFile,'\n');
        end
        
        fclose(EdgesMatrixFile)
    
        % Display the size of the resulting table
        disp('Valid edges found:');
        size(valid_edges)
    else
        disp('Failed to find valid edges within the maximum number of attempts.');
        return
    end

end

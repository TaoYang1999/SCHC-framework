function regionStats = analyzeSurface(surfaceMatrix, Z)
% Input:
% Surfacematrix - matrix representing the height of rough surface topography
% Z - height threshold, below which the element is Nan
% Output:
% Regionstats - cell array, containing statistics of each region

global delta_L

surfaceMatrix(surfaceMatrix < Z) = NaN;

binaryMatrix = ~isnan(surfaceMatrix); 

cc = bwconncomp(binaryMatrix); 

stats = regionprops(cc, 'Area', 'PixelIdxList', 'BoundingBox', 'Centroid', 'Perimeter');

numRegions = numel(stats);
regionStats = cell(numRegions, 6); 

boundaries = bwboundaries(binaryMatrix);  

for i = 1:numRegions
    %% Column 1: area Index
    regionStats{i, 1} = i;

    %% Column 2: calculated area
    regionStats{i, 2} = stats(i).Area * delta_L^2; 

    %% Column 3: get the pixel index of the current region and calculate the maximum height
    heightValues = surfaceMatrix(stats(i).PixelIdxList);
    regionStats{i, 3} = max(heightValues);

    %% Column 4: get centroid coordinates
    centroid = stats(i).Centroid; 
    regionStats{i, 4} = centroid * delta_L;

    %% Column 5: perimeter of area.    NOTE: only outer boundary
    regionStats{i, 5} = stats(i).Perimeter * delta_L;  

    %% Column 6: circularity of area    1
    boundaryPoints = boundaries{i} * delta_L;
    boundaryPoints = [boundaryPoints(:,2), boundaryPoints(:,1)];
    centroidX = centroid(1) * delta_L; 
    centroidY = centroid(2) * delta_L;  
    distances = sqrt((boundaryPoints(:,1) - centroidX).^2 + (boundaryPoints(:,2) - centroidY).^2);
    regionStats{i, 6} = regionStats{i, 2}/(max(distances)^2*pi);  

    %% Column 7: circularity of area   2
    boundaries = bwboundaries(binaryMatrix, 'noholes');  
    regionStats{i, 7} = 4*pi*regionStats{i, 2}/regionStats{i, 5}^2;    % CIR = (4 * pi * A) / (P^2);
end

end

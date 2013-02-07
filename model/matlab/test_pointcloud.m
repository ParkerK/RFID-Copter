% close previous figures
close all
clear all
pause on
hold on

% constants
XMIN = -20;
XMAX = 20;
YMAX = 20;
ZMIN = -2;
ZSTEP = 0.5;
ZMAX = 2;
NUMZSTEPS = (ZMAX - ZMIN) / ZSTEP + 1;
POINT_SIZE = 10;
PTS_PER_SLICE = 100;

model_3d = [];

% create model
for k = 0:(NUMZSTEPS-1)
    
    % x y z intensity
    slice_2d = zeros(PTS_PER_SLICE,4);
    
    % z coordinate
    slice_2d(:,3) = ZMIN + ZSTEP * k;
    
    % x coordinate
    slice_2d(:,1) = XMIN + rand(PTS_PER_SLICE,1) * 2 * XMAX;
    
    % y coordinate
    slice_2d(:,2) = 0 + rand(PTS_PER_SLICE,1) * YMAX;
    
    % linear intensity model
    dist = sqrt( slice_2d(:,1) .^ 2 + slice_2d(:,2) .^ 2);
    max_dist = sqrt( XMAX ^ 2 + YMAX ^ 2);
    slice_2d(:,4) = (1 - dist ./ max_dist);
    
    % concatenate to 3d matrix for this update step
    model_3d = vertcat(model_3d, slice_2d);
    
end

% plot it
set(gcf,'OuterPosition',[1120 0 800 800])
scatter3(model_3d(:,1), model_3d(:,2), model_3d(:,3), POINT_SIZE, model_3d(:,4));

xlabel('X (meters)')
ylabel('Y (meters)')
zlabel('Z (meters)')
campos([-10, -5, 5])
% close previous figures
close all
clear all
clc
pause on
hold on

%% Create Model
% constants
XMIN = -10;
XMAX = 10;
YMIN = 0;
YMAX = 20;
ZMIN = 0;
ZSTEP = 0;
ZMAX = 0;
NUMZSTEPS = (ZMAX - ZMIN) / ZSTEP + 1;
POINT_SIZE = 10;

% model readings
% model2d [ x y prob ]
pts_per_row = (XMAX-XMIN)*2+1;
total_pts = pts_per_row * (YMAX + 1);
model_2d = [];
max_dist = sqrt( XMAX ^ 2 + YMAX ^ 2);
for y = 0:YMAX
    for x = XMIN:XMAX
        prob = 1 - sqrt( x ^ 2 + y ^ 2) / max_dist;
        model_2d = vertcat(model_2d, [x y prob]);
    end
end

set(gcf, 'Position', [1300 0 800 800])
scatter( model_2d(:,1), model_2d(:,2), 10, model_2d(:,3) )
xlabel('X (meters)')
ylabel('Y (meters)')

%% Particle Filter
NUM_PARTICLES = 1000;
NOISE_SIGMA_MEASUREMENT_PROB = 0.2;
NOISE_SIGMA_RESAMPLE = 0.5;

% model input data
target = 0.5;
num_readings = 10;
readings = zeros(num_readings,1);
for k = 1:num_readings
    readings(k,1) = target;
end

% randomly initialize them
% w [ x y error norm]
w = zeros(NUM_PARTICLES,4);
% randomize intially, not gaussian because this needs to be random over
% everything
w(:,1) = XMIN + rand(NUM_PARTICLES,1) * 2 * XMAX;
w(:,2) = 0 + rand(NUM_PARTICLES,1) * YMAX;
f = scatter( w(:,1), w(:,2), 'k', 'filled');

for k = 1:size(readings,1);
    
    % find probability measurement came from a given particle
    for i = 1:NUM_PARTICLES
        % find RSSI value that corresponds to w
        % match wx,wy with the RSSI value at model2dx, model2dy
        % round it to then nearest one
        wx = round(w(i,1));
        wy = round(w(i,2));
                
        model_index = intersect( find(model_2d(:,1)== wx), find(model_2d(:,2)==wy) );
        
        
        w(i,3) = normpdf( abs( model_2d(model_index,3) - readings(k) ), 0, NOISE_SIGMA_MEASUREMENT_PROB );
%         
%         model_value = model_2d(model_index,3)
%         abs_val_diff = abs( model_2d(model_index,3) - readings(k) )       
%         prob_measurement = w(i,3)
        
        
    end
        
    sum_of_errors = sum(w(:,3))
        
    for i = 1:NUM_PARTICLES
        w(i,4) = w(i,3) / sum_of_errors;
    end
    
    % Resampling
    w2 = zeros(NUM_PARTICLES,4);
    % initialize index
    index = randi(NUM_PARTICLES);
    beta = 0.0;
    mw = max(w(:,4));
    for i = 1:NUM_PARTICLES
        beta = beta + rand * 2.0 * mw;
        while beta > w(index,4)
            beta = beta - w(index,4);
            index = mod(index + 1, NUM_PARTICLES)+1;
        end
        
        w2(i,1) = w(index,1) + normrnd(0,NOISE_SIGMA_RESAMPLE);
        w2(i,2) = w(index,2) + normrnd(0,NOISE_SIGMA_RESAMPLE);
        w2(i,3) = w(index,3);
        w2(i,4) = w(index,4);
        
        if(w2(i,1) > XMAX)
            w2(i,1) = XMAX;
        end
        if(w2(i,1) < XMIN)
            w2(i,1) = XMIN;
        end
        if(w2(i,2) > YMAX)
            w2(i,2) = YMAX;
        end
        if(w2(i,2) < YMIN)
            w2(i,2) = YMIN;
        end
    end

    w = w2;
 
    pause(0.3)
    delete(f)
    f = scatter( w(:,1), w(:,2), 'k', 'filled' );

    
end




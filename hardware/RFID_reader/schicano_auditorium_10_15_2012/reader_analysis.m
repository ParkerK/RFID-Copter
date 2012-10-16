clear all
close all

FigHandle = figure;
%set(FigHandle, 'Position', [800, 200, 1049, 895]);

%% AVERAGE RSSI PLOT

%subplot(2,1,1)

% this is the matrix of averaged rssi values across 2D space
avg_rssi_matrix_right = zeros(30,20);

for i = 1:26
    % import data from files
    filename = strcat('logs/',int2str(i),'ft_c_to_l.txt');
    data = importdata(filename);
    
    % replace NaN's with 0's
    data(isnan(data)) = 0;

    % number of rows of data in each file varies
    [rows,cols] = size(data);
    for j = 1:rows
        % 30 readings at each point, take average of that
        avg_rssi = mean(data(j,2:length(data)));
        avg_rssi_matrix_right(i,j) = avg_rssi;
    end       
end

% make left matrix
% data was actually taken on the left
avg_rssi_matrix_left = fliplr(avg_rssi_matrix_right);
% concatenate matrices
avg_rssi_matrix = horzcat(avg_rssi_matrix_left, avg_rssi_matrix_right);
% plot 
filtered_data = imfilter(avg_rssi_matrix, fspecial('average',[2 2]));
avg_rssi_plot = surf(filtered_data)
title('Average RSSI Moving Over 2D Space  ', 'FontSize', 40)
zlabel('Average RSSI  ','FontSize',20)
xlabel('X in feet  ','FontSize',20)
ylabel('Y in feet  ','FontSize',20)

%% PERCENT OF READER MISSES PLOT

figure
%subplot(2,1,2)

% this is the matrix of averaged rssi values across 2D space
x = 30;
y = 20
reader_misses_matrix_right = zeros(x,y);
for i = 1:x
    for j = 1:y
        reader_misses_matrix_right(i,j) = 100;
    end
end

for i = 1:26
    % import data from files
    filename = strcat('logs/',int2str(i),'ft_c_to_l.txt');
    data = importdata(filename);
    
    % replace NaN's with 0's
    data(isnan(data)) = 0;
    
    % number of rows of data in each file varies
    [rows,cols] = size(data);
    for j = 1:rows
        % 30 readings at each point, find how many were misses
        reader_misses = sum(data(j,2:length(data)) == 0);
        reader_misses_matrix_right(i,j) = reader_misses / 30 * 100;
    end  
end

% make left matrix
% data was actually taken on the left
reader_misses_matrix_left = fliplr(reader_misses_matrix_right);
% concatenate matrices
reader_misses_matrix = horzcat(reader_misses_matrix_left, reader_misses_matrix_right);
% plot 
filtered_data = imfilter(reader_misses_matrix, fspecial('average',[2 2]));
% fix edges from the filter
for i = 1:30
    filtered_data(i,40) = 100;
end
for i = 1:40
    filtered_data(30,i) = 100;
end   
% plot 
reader_misses_plot = surf(filtered_data)
title('Percent of Reader Misses Over 2D Space  ','FontSize',40)
zlabel('Percent of Misses  ','FontSize',20)
xlabel('X in feet  ','FontSize',20)
ylabel('Y in feet  ','FontSize',20)
colorbar



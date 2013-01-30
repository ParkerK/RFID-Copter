clear all
close all

FigHandle = figure;
set(FigHandle, 'Position', [800, 200, 1049, 895]);

titles = {'Average RSSI at 50.5 Inches - Level With Reader'; 'Average RSSI at 70 Inches'; 'Average RSSI at 89 Inches'; 'Average RSSI at 130 Inches'};
filenames = {'v50.5_d';'v70_d';'v89_d';'v130_d'};

for k = 1:4
    
    subplot(2,2,k)

    % this will be the matrix of averaged rssi values across 2D space
    %avg_rssi_matrix_right = zeros(10,7);
    avg_rssi_matrix_right = zeros(31,46);

    % we took data every 5 feet out to 45 feet, thats 10 total places
    for i = 0:9
        % import data from files
        filename = strcat('logs/',filenames{k},int2str(i*5),'ft.txt');
        data = importdata(filename);

        % get num of rows / columns
        [rows,cols] = size(data);
        
        for j = 1:rows
            % 30 readings at each point, take average of that
            avg_rssi = mean(data(j,2:cols));
            avg_rssi_matrix_right(1+((j-1)*5),1+(i*5)) = avg_rssi
            
        end       
    end

    avg_rssi_matrix = zeros(62,46);
    avg_rssi_matrix = [zeros(31,46); avg_rssi_matrix_right];
    avg_rssi_matrix = [fliplr(fliplr(flipud(avg_rssi_matrix_right))); avg_rssi_matrix_right];

    % plot 
    %filtered_data = imfilter(avg_rssi_matrix, fspecial('average',[2 2]));
    rssi_plot = surf(avg_rssi_matrix)
    title(titles{k})
    zlabel('Average RSSI')
    xlabel('X in feet')
    ylabel('Y in feet')
    
    colorbar
end


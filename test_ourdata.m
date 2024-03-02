%% Version: 02/Feb/2024, Writen and Edited by Yanan
   % Use the Pan Tompking Algorithm to find Peaks Point R
   % Based on the obtained Peaks Point R to find the S Point
   clear;
   clc;
   set(0,'defaultfigurecolor','w') 

   % Load the data from CSV file
data = readtable('AFE4500_CAPTURED_DATA_copy.csv');

% Extract the columns by their names
t_ecg = data.ECG_Time/4;
ecg = data.ECG_Value;
t_ppg = data.TIA1_3_Time*2.5;
ppg = data.TIA1_3_Value;

nonNaNIndex = ~isnan(ppg);
ppg = ppg(nonNaNIndex);
t_ppg = t_ppg(nonNaNIndex);
%% Something need to be edited is:
%  1. the search_window of S point -- 129
%  2. The find_peaks function in pan_tompkin.m to find the peak point of -- 92 
%  3. window_size for detecting the Lowest point of PPG wave -- 225
%  4. find_peaks for PPG peak edit the MinPeakDistance --213
%  5. Edit the chosen part of ECG and PPG for extract features -- 360-361
%%
% %% Load Captured Signal
% % Load data from CSV file
figure
plot(t_ecg,ecg)
figure
plot(t_ppg,ppg)
% %% Load data from mat file
% % 加载.mat文件，这里假设data结构体已经在工作空间中
% data = load('part_5.mat');
% % 确定cell数组的长度
% numCells = length(data.p);
% 
% % 遍历cell数组
% for i = 1:numCells
%     % 从data结构体中提取第i个cell的内容
%     cellContent = data.p{i};
% 
%     % 构建文件名，包含cell的索引
%     filename = sprintf('p_Record_%d.mat', i);
% 
%     % 保存提取的数据到.mat文件
%     save(filename, 'cellContent');
% end
% 
% % 假设我们关注的是data.p中的第一个cell
% cellContent = data.p{6}; % 选择第一个cell的内容
% 
% % 获取采样率
% Fs = 125; % 采样率为125Hz
% t = (0:length(cellContent)-1) / Fs; % 创建时间向量
% 
% % 绘制第1行数据（例如，PPG信号）
% subplot(3, 1, 1); % 创建3行1列的子图，这是第1个
% plot(t, cellContent(1, :));
% title('PPG Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% 
% % 绘制第2行数据（例如，ABP信号）
% subplot(3, 1, 2); % 第2个子图
% plot(t, cellContent(2, :));
% title('ABP Signal');
% xlabel('Time (s)');
% ylabel('mmHg');
% 
% % 绘制第3行数据（例如，ECG信号）
% subplot(3, 1, 3); % 第3个子图
% plot(t, cellContent(3, :));
% title('ECG Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');


%% Apply Savitzky-Golay Filter to Captured Signal
% fs = 125;
% ecg = cellContent(3, :);
% t_ecg = (0:length(ecg) - 1)/fs;
% ppg = cellContent(1, :);
% t_ppg = (0:length(ppg) - 1)/fs;

% windowSize_ECG = 51;  % Window size for the filter (should be odd)
% polyOrder_ECG = 3;   % Polynomial order
% 
% % Apply Savitzky-Golay Filter to ECG and PPG data
% filtered_ecg = sgolayfilt(ecg, polyOrder_ECG, windowSize_ECG);


%% Do Pan Tompking Algortihm to find Peaks Point R
fs = 500;
time_difference = 0.8;
 [qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecg,fs,1,time_difference);
time_diff = diff(qrs_i_raw)/(fs);

% Calculate the second derivative of the ECG signal
second_derivative_ecg = diff(ecg, 2);

% Initialize array for refined R points
% refined_qrs_i = zeros(size(qrs_i_raw));
% 
% for i = 1:length(qrs_i_raw)
%     current_idx = qrs_i_raw(i);  % Start from the initial peak index
%     peak_found = false;
% 
%     while ~peak_found && current_idx < (length(ecg) - 1)
%         % Condition to check if the current point is a peak
%         is_peak = (ecg(current_idx) > ecg(current_idx - 1)) && (ecg(current_idx) > ecg(current_idx + 1));
% 
%         % Ensure we're not at the boundaries where second_derivative_ecg might be out of bounds
%         if current_idx <= length(second_derivative_ecg) + 1
%             % Additional condition: current point must be a local maximum (second_derivative < 0)
%             is_peak = is_peak && (second_derivative_ecg(current_idx - 1) < 0);
%         end
% 
%         if is_peak
%             peak_found = true;
%             refined_qrs_i(i) = current_idx;
%         else
%             % Move to the next point
%             current_idx = current_idx + 1;
%         end
%     end
% end

% Initialize array for S points
s_points = zeros(size(qrs_i_raw));

% Define search window for S point (example: 0.5 seconds after R peak)
search_window = 1.2; % 0.5 s, adjust as necessary
num_samples_in_window = round(search_window * fs);

% Loop through R peaks to find S points
for i = 1:length(qrs_i_raw)
    % Define the search window for the current R peak
    start_idx = qrs_i_raw(i);
    end_idx = min(start_idx + num_samples_in_window, length(ecg));
    
    % Search for local minima in the window
    for j = start_idx+1:end_idx-1
        if (ecg(j) < ecg(j-1)) && (ecg(j) < ecg(j+1))
            if second_derivative_ecg(j-1) > 0
                s_points(i) = j;
                break;
            end
        end
    end
end


% % Define the Index Region
% startIndex = 1;
% endIndex = 104;
% 
% endIndex = min(endIndex, length(qrs_i_raw));

% Select 100 S and R points

selected_R_points = qrs_i_raw;
selected_S_points = s_points;

% Select Corresponding time to plot
selected_R_times = t_ecg(selected_R_points);
selected_R_amplitudes = ecg(selected_R_points);
selected_S_times = t_ecg(selected_S_points);
selected_S_amplitudes = ecg(selected_S_points);

% Plot the selected R S points
figure;
plot(t_ecg, ecg);
hold on;
plot(selected_R_times, selected_R_amplitudes, 'or', 'MarkerSize', 8); % 选中的R点
plot(selected_S_times, selected_S_amplitudes, 'og', 'MarkerSize', 8); % 选中的S点
hold off;
title('Filtered ECG with Selected R and S Points');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Filtered ECG', 'Selected R Points', 'Selected S Points');
% Obtain the time difference
tdiff_R = diff(selected_R_times);
tdiff_S =diff(selected_S_times);
tdiff_SR = selected_S_times - selected_R_times;
%% Print the final results
fprintf("The mean of time difference of R peaks(ECG) is: %f\n", mean(tdiff_R));
fprintf("The variance of time difference of R peaks(ECG) is: %f\n", var(tdiff_R));
fprintf("The mean of time difference of S peaks(ECG) is: %f\n", mean(tdiff_S));
fprintf("The variance of time difference of S peaks(ECG) is: %f\n", var(tdiff_S));
fprintf("The mean of time difference of SR peaks(ECG) is: %f\n", mean(tdiff_SR));
fprintf("The variance of time difference of SR peaks(ECG) is: %f\n", var(tdiff_SR));


%% Load PPG data from CSV file
%Flip Signal
ppg = -ppg;

% Detrend the signal
ppg_detrended = detrend(ppg);

% Bandpass filter design
Fs = 50; % Replace with your actual sampling frequency
Wp = [0.5 4] / (Fs/2); % Passband Frequencies (normalized)
Ws = [0.1 10] / (Fs/2); % Stopband Frequencies (normalized)
[Rp, Rs] = deal(3, 30); % Passband ripple and Stopband attenuation
[n, Wn] = buttord(Wp, Ws, Rp, Rs); % Filter order and cutoff frequency
[b, a] = butter(n, Wn, 'bandpass'); % Filter coefficients

% Apply bandpass filter
ppg_filtered = filtfilt(b, a, ppg_detrended);

% Optional: Smooth the data
ppg_smooth = smoothdata(ppg_filtered, 'movmean', 10); % Window size of 10, adjust as needed

% Find peaks
[peaks, locs] = findpeaks(ppg_smooth, 'MinPeakHeight', mean(ppg_smooth), 'MinPeakDistance', 20); % Adjust threshold and distance

% Plot the results
figure;
plot(t_ppg, ppg_smooth);
hold on;
plot(t_ppg(locs), peaks, 'ro');
title('Filtered PPG Signal with Peaks');
xlabel('Time (s)');
ylabel('Amplitude');
hold off;
% Define a window size for local minimum search
window_size = round(40); % Adjust this value as needed

% Calculate the second derivative
second_derivative = diff(ppg_smooth, 2);

% Initialize vector to hold the indices of the minima
minima_indices = zeros(size(locs));

% Find minima before peaks with windowed search
for i = 1:length(locs)
    % Ensure the window does not go beyond the start of the signal
    search_start = max(locs(i) - window_size, 1);
    search_end = locs(i) - 1; % Stop search just before the peak
    
    % Get the local window to search for the minimum
    local_window = ppg_smooth(search_start:search_end);
    
    % Find the index of the minimum value in this window
    [~, min_idx] = min(local_window);
    
    % Correct the index relative to the entire signal
    minima_indices(i) = search_start + min_idx - 1;
    
    % Ensure that the second derivative at this point is greater than zero
    % if second_derivative(minima_indices(i) - 1) <= 0
    %     minima_indices(i) = []; % If condition not met, mark as NaN
    % end
end

% Clean up NaN entries from the minima indices
minima_indices = minima_indices(~isnan(minima_indices));

% Plot the results with minima points
figure;
plot(t_ppg, ppg_smooth);
hold on;
plot(t_ppg(locs), peaks, 'ro', 'MarkerFaceColor', 'r'); % Peaks marked with red circles
for i = 1:length(minima_indices)
plot(t_ppg(minima_indices(i)), ppg_smooth(minima_indices(i)), 'ko', 'MarkerFaceColor', 'k'); % Minima marked with black circles
end
title('Filtered PPG Signal with Peaks and Minima');
xlabel('Time (s)');
ylabel('Amplitude');
hold off;

%% Define the Index Region to find selected points of PPG
% startIndex = 5;
% endIndex = 104;


% Select 100 S and R points
selected_R_points_PPG = locs;
selected_S_points_PPG = minima_indices;
% Obtain the time difference
tdiff_R_PPG = diff(selected_R_points_PPG)/Fs;
tdiff_S_PPG =diff(selected_S_points_PPG)/Fs;
tdiff_SR_PPG = (selected_R_points_PPG - selected_S_points_PPG)/Fs;

% Plot the results with selected points
figure;
plot(t_ppg, ppg_smooth);
hold on;
plot(t_ppg(selected_R_points_PPG), ppg_smooth(selected_R_points_PPG), 'ro', 'MarkerFaceColor', 'r'); % Peaks marked with red circles
for i = 1:length(selected_S_points_PPG)
plot(t_ppg(selected_S_points_PPG(i)), ppg_smooth(selected_S_points_PPG(i)), 'ko', 'MarkerFaceColor', 'k'); % Minima marked with black circles
end
title('Selected Peaks and Minima Points of PPG signals');
xlabel('Time (s)');
ylabel('Amplitude');
hold off;

% Print the time difference characteristics of Peaks and minimum points of PPG signal
fprintf("The mean of time difference of R peaks(PPG) is: %f seconds \n", mean(tdiff_R_PPG));
fprintf("The variance of time difference of R peaks(PPG) is: %f seconds \n", var(tdiff_R_PPG));
fprintf("The mean of time difference of S peaks(PPG) is: %f seconds \n", mean(tdiff_S_PPG));
fprintf("The variance of time difference of S peaks(PPG) is: %f seconds \n", var(tdiff_S_PPG));
fprintf("The mean of time difference of SR peaks(PPG) is: %f seconds \n", mean(tdiff_SR_PPG));
fprintf("The variance of time difference of SR peaks(PPG) is: %f seconds \n", var(tdiff_SR_PPG));

%% Find time difference between PPG peaks and ECG peaks
% shift the ECG's peak points
% startIndex = 5;
% endIndex = 104;
% 
% endIndex = min(endIndex, length(qrs_i_raw));

% Select 100 S and R points
selected_R_points = qrs_i_raw;

% Initialize a variable to store the time differences
time_diff_ECG_PPG = zeros(length(selected_R_points), 1);

% Loop over each R peak to find the second PPG peak following it
for i = 1:length(selected_R_points)
    R_peak_time = t_ecg(selected_R_points(i));
    
    % Find all PPG peaks after the R peak
    subsequent_PPG_peak_indices = find(t_ppg(locs) > R_peak_time);
    
    % Check if there are at least two PPG peaks after the R peak
    if length(subsequent_PPG_peak_indices) >= 1
        % If so, select the second PPG peak
        second_PPG_peak_index = subsequent_PPG_peak_indices(1);
        % Get the time of the second PPG peak
        PPG_peak_time = t_ppg(locs(second_PPG_peak_index));
        % Calculate the time difference between the R peak and the second PPG peak
    time_diff_ECG_PPG(i) = PPG_peak_time - R_peak_time;
    else
    % If there are less than two PPG peaks after the R peak, set to NaN
    time_diff_ECG_PPG(i) = NaN;
    end
end

% Calculate the mean and variance of the time differences
mean_time_diff = mean(time_diff_ECG_PPG);
var_time_diff = var(time_diff_ECG_PPG);

% Display the results
fprintf('The mean time difference between ECG R-peaks and subsequent PPG peaks is %.3f seconds.\n', mean_time_diff);
fprintf('The variance of time differences between ECG R-peaks and subsequent PPG peaks is %.5f seconds squared.\n', var_time_diff);

%% Find time difference between PPG lowest points and ECG peaks
% Loop over each R peak to find the next lowest PPG point
% shift the ECG's peak points
% startIndex = 4;
% endIndex = 103;
% 
% endIndex = min(endIndex, length(qrs_i_raw));

% Select 100 S and R points
selected_R_points = qrs_i_raw;
selected_S_points = s_points;

% Select Corresponding time to plot (2 data difference)
selected_R_times = t_ecg(selected_R_points(10:50));
selected_ppg_times = t_ppg(selected_S_points_PPG(11:51));

lowest_point_time_diffs = selected_ppg_times - selected_R_times;

% Display the results
fprintf('The mean time difference between ECG R-peaks and the next lowest PPG points is %.3f seconds.\n', mean(lowest_point_time_diffs));
fprintf('The variance of time differences between ECG R-peaks and the next lowest PPG points is %.5f seconds squared.\n', var(lowest_point_time_diffs));

result = [mean(tdiff_R),mean(tdiff_S),mean(tdiff_SR),mean(tdiff_R_PPG),mean(tdiff_SR_PPG),mean_time_diff,mean(lowest_point_time_diffs);var(tdiff_R),var(tdiff_S),var(tdiff_SR),var(tdiff_R_PPG),var(tdiff_SR_PPG),var_time_diff,var(lowest_point_time_diffs)];
RowName = {'Mean';'Variance'};
ColumnNames = {'RR_ECG','SS_ECG','RS_ECG','RR_PPG','RL_PPG','RECG_RPPG','RECG_LPPG'};
savepath = 'C:\Users\Administrator\Desktop\ISL Lab\Lab_Data\Matlab_data\Test_1_26\after-eating-chips-20240126T212523Z-001\after-eating-chips\result.xlsx';
kWriteTable(result,RowName,ColumnNames,savepath);
%%
% 确定最大导数点
ppg_diff1 = diff(ppg_smooth);
[max_deriv_values, max_deriv_locs] = findpeaks(ppg_diff1,'MinPeakDistance', 40); % 寻找PPG信号导数的峰值点
% 初始化存储IT点的数组
it_points_x = [];
it_points_y = [];

% 遍历所有谷点
for i = 15:30
    % 当前谷点的位置和导数（斜率）
    min_loc = minima_indices(i);
    min_slope = ppg_diff1(min_loc); % 谷点处切线斜率可以近似为0
    min_intercept = ppg_smooth(min_loc); % y截距等于谷点的PPG值
    
    % 找到谷点后的第一个最大导数点
    max_deriv_loc = max_deriv_locs(find(max_deriv_locs > min_loc, 1, 'first'));
    
    if isempty(max_deriv_loc)
        continue; % 如果找不到最大导数点，则跳过当前谷点
    end
    
    % 最大导数点的位置和导数（斜率）
    max_deriv_slope = ppg_diff1(max_deriv_loc);
    max_deriv_intercept = ppg_smooth(max_deriv_loc) - max_deriv_slope * t_ppg(max_deriv_loc);
    
    % 计算交点
    % 两条切线方程为 y = min_slope * x + min_intercept 和 y = max_deriv_slope * x + max_deriv_intercept
    % 解这两个方程找到交点 x 坐标
    it_x = (min_intercept - max_deriv_intercept) / (max_deriv_slope - min_slope);
    it_y = max_deriv_slope * it_x + max_deriv_intercept;
    
    % 存储交点
    it_points_x(end + 1) = it_x;
    it_points_y(end + 1) = it_y;
end

figure;
plot(t_ppg, ppg_smooth, 'b'); % 绘制PPG信号
hold on;

% 标记最小值点
plot(t_ppg(minima_indices), ppg_smooth(minima_indices), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 

% 标记最大导数点
plot(t_ppg(max_deriv_locs), ppg_smooth(max_deriv_locs), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8); 

% 标记交叉点（IT点）
plot(it_points_x, it_points_y, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); 

hold off;
title('PPG Signal with IT Points');
xlabel('Time (s)');
ylabel('Amplitude');
legend('PPG Signal', 'Valley Points', 'Max Derivative Points', 'IT Points');


for i = 1:length(it_points_x)
    t_ppg_lowest = t_ppg(minima_indices);
    % Find all ppg minima before it_points
    ppg_minima_before_it = find(it_points_x(i)>t_ppg_lowest);
    PAT_PPG4 = it_points_x - t_ppg_lowest(end);
end
fprintf("The mean of PAT_PPG4 is %.3f.\n",mean(PAT_PPG4));
fprintf("The variance of PAT_PPG4 is %.5f.\n",var(PAT_PPG4));
%%
figure;
hold on;

% 绘制ECG信号
plot(t_ecg, ecg, 'b', 'LineWidth', 1); % 使用蓝色绘制ECG
% 绘制PPG信号，可能需要调整PPG信号的幅度以便在同一图中更好地显示
plot(t_ppg, ppg_smooth * max(abs(ecg)) / max(abs(ppg_smooth)), 'm', 'LineWidth', 1); % 使用品红色绘制PPG

% 标注ECG的R点
selected_R_times = t_ecg(selected_R_points);
selected_R_amplitudes = ecg(selected_R_points);
plot(selected_R_times, selected_R_amplitudes, 'or', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % 使用红色圆圈标注

% 标注PPG的R点
plot(t_ppg(selected_R_points_PPG), ppg_smooth(selected_R_points_PPG) * max(abs(ecg)) / max(abs(ppg_smooth)), 'xr', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % 使用红色叉号标注

% 标注ECG的S点
selected_S_times = t_ecg(selected_S_points);
selected_S_amplitudes = ecg(selected_S_points);
plot(selected_S_times, selected_S_amplitudes, 'og', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % 使用绿色圆圈标注

% 标注PPG的S点
plot(t_ppg(selected_S_points_PPG), ppg_smooth(selected_S_points_PPG) * max(abs(ecg)) / max(abs(ppg_smooth)), 'xg', 'MarkerSize', 8, 'MarkerFaceColor', 'g'); % 使用绿色叉号标注

hold off;

title('ECG and PPG Signals with Selected R and S Points');
xlabel('Time (s)');
ylabel('Amplitude');
legend('ECG', 'PPG', 'ECG R Points', 'PPG R Points', 'ECG S Points', 'PPG S Points');
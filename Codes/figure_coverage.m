clear all;
close all;

% Aggregate the Monte Carlo simulations of the coverage simulations and
% create the summary figures comparing the CBF case with the baseline.

% 
current_path = pwd;
path = fullfile(current_path, 'heat-ef=001_5o');
for i = 1:50
    % Load one realization of the safety-filtered dynamics.
    filename = strcat("data", num2str(i), '.mat');
    path_file = fullfile(path, num2str(i), filename);
    data = open(path_file);
    nf_inside_Rk(i,:) = data.nf_inside;
    Hf(i,:) = data.Hf;
end

% Ensemble average of the obstacle-incursion fraction under CBF.
nf_inside_Rk_mean = mean(nf_inside_Rk,1) / 720;

path = fullfile(current_path, 'heat-nocbf_5o');
for i = 1:50
    % Load one realization of the baseline dynamics without safety filter.
    filename = strcat("data", num2str(i), '.mat');
    path_file = fullfile(path, num2str(i), filename);
    data = open(path_file);
    nf_inside_nocbf(i,:) = data.nf_inside;
    Hf_nocbf(i,:) = data.Hf;
end



% Ensemble average of the obstacle-incursion fraction without CBF.
nf_inside_nocbf_mean = mean(nf_inside_nocbf,1) / 720;


t = 0:0.01:50;
column_width = 7*0.49*0.75;   % Column width in inches
aspect_ratio = 0.6;     % height / width
figure_height = column_width * aspect_ratio;

nf_inside_Rk_mean2 = sum(nf_inside_Rk,2)*0.01;
nf_inside_nocbf_mean2 = sum(nf_inside_nocbf,2)*0.01;

% Boxplot of the time-integrated obstacle occupancy.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
set(gca,'FontSize',10,'FontName',"Times New Roman")
set(gca,'LooseInset',get(gca,'TightInset'))
axis tight
boxchart(2*ones(size(nf_inside_Rk_mean2)), nf_inside_Rk_mean2);          % x = 1
boxchart(ones(size(nf_inside_nocbf_mean2)), nf_inside_nocbf_mean2);        % x = 2
xticks([1 2])
xticklabels({'NO CBF', 'CBF'})
ylabel('$\int_T N^{\mathrm{in}} \, dt$', 'Interpreter', 'latex')

nf_inside_Rk_std = std(nf_inside_Rk, 0, 1) / 720;
nf_inside_nocbf_std = std(nf_inside_nocbf, 0, 1) / 720;



%% Follower
% Mean +/- one standard deviation of the obstacle-incursion fraction.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
% Ensure proper scaling for LaTeX documents
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
set(gca,'FontSize',10,'FontName',"Times New Roman")
set(gca,'LooseInset',get(gca,'TightInset'))
axis tight
% Configure tick labels with LaTeX interpreter and fontsize
h1 = fill([t fliplr(t)], ...
     [nf_inside_nocbf_mean + nf_inside_nocbf_std, fliplr(nf_inside_nocbf_mean - nf_inside_nocbf_std)], ...
     [0.0660 0.4430 0.7450], 'EdgeColor', [0.0660 0.4430 0.7450], 'FaceAlpha', 0.3);
plot(t, nf_inside_nocbf_mean, 'Color', '[0.0660 0.4430 0.7450]', 'LineWidth', 1.2);
h2 = fill([t fliplr(t)], ...
     [nf_inside_Rk_mean + nf_inside_Rk_std, fliplr(nf_inside_Rk_mean - nf_inside_Rk_std)], ...
     [0.8510 0.3255 0.0980], 'EdgeColor', [0.8510 0.3255 0.0980], 'FaceAlpha', 0.3);
plot(t, nf_inside_Rk_mean, 'Color', [0.8510 0.3255 0.0980], 'LineWidth', 1.2);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\langle N^{\mathrm{in}}\rangle / N$', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');
lgd = legend([h1, h2], {'NO CBF', 'CBF'},  'Orientation','horizontal', ...
    'Location','southeast');
lgd.FontSize = 8;
lgd.ItemTokenSize = [8 6];
%% Export all figures to PDF

% Export every open figure as a vector PDF.
fig_folder = fullfile(pwd,'figure');

if ~exist(fig_folder,'dir')
    mkdir(fig_folder);
end

figs = findall(groot,'Type','figure');

for k = 1:length(figs)

    fig = figs(k);

    filename = fullfile(fig_folder, sprintf('figure_heat_%02d.pdf',k));

    exportgraphics(fig, filename, ...
        'ContentType','vector', ...
        'BackgroundColor','none');

end

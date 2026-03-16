clear all;
close all;

% Aggregate the Monte Carlo simulations of the leader-follower experiment
% and generate the summary figures used in the paper.

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultTextInterpreter','latex');          
set(groot,'defaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize', 16);
set(groot, 'DefaultTextFontSize', 16);
set(groot, "defaultFigurePosition", [0 0 560 420])
set(groot, "defaultFigureWindowStyle", "normal")


% 
current_path = pwd;
path = fullfile(current_path, 'el=0013-ef=001-5o');
for i = 1:50
    % Load one realization and collect the time histories needed for statistics.
    filename = strcat("data", num2str(i), '.mat');
    path_file = fullfile(path, num2str(i), filename);
    data = open(path_file);
    kl_leader_Rk(i,:) = data.dkl_l_hist;
    kl_follower_Rk(i,:) = data.dkl_f_hist;
    nl_inside_Rk(i,:) = data.nl_inside;
    nf_inside_Rk(i,:) = data.nf_inside;
    Hl(i,:) = data.Hl;
    chi(i,:) = data.chi/720;
    Hf(i,:) = data.Hf;
end

% Compute ensemble averages for the CBF-controlled case.
kl_leader_Rk_mean = mean(kl_leader_Rk,1);
kl_follower_Rk_mean = mean(kl_follower_Rk,1);
nl_inside_Rk_mean = mean(nl_inside_Rk,1) / 100;
nf_inside_Rk_mean = mean(nf_inside_Rk,1) / 720;

path = fullfile(current_path, 'leader-follower_nocbf_5o');
for i = 1:50
    % Load the corresponding realization of the baseline without CBF.
    filename = strcat("data", num2str(i), '.mat');
    path_file = fullfile(path, num2str(i), filename);
    data = open(path_file);
    kl_leader_nocbf(i,:) = data.dkl_l_hist;
    kl_follower_nocbf(i,:) = data.dkl_f_hist;
    nl_inside_nocbf(i,:) = data.nl_inside;
    nf_inside_nocbf(i,:) = data.nf_inside;
    Hl_nocbf(i,:) = data.Hl;
    chi_nocbf(i,:) = data.chi/720;
    Hf_nocbf(i,:) = data.Hf;
end

% Compute ensemble averages for the uncontrolled/baseline case.
kl_leader_nocbf_mean = mean(kl_leader_nocbf,1);
kl_follower_nocbf_mean = mean(kl_follower_nocbf,1);
nl_inside_nocbf_mean = mean(nl_inside_nocbf,1) / 100;
nf_inside_nocbf_mean = mean(nf_inside_nocbf,1) / 720;


t = 0:0.01:100;
column_width = 7*0.49*0.49; % Column width in inches
aspect_ratio = 0.75;     % height / width
figure_height = column_width * aspect_ratio;

nl_inside_Rk_mean2 = sum(nl_inside_Rk,2)*0.01;
nf_inside_Rk_mean2 = sum(nf_inside_Rk,2)*0.01;
nl_inside_nocbf_mean2 = sum(nl_inside_nocbf,2)*0.01;
nf_inside_nocbf_mean2 = sum(nf_inside_nocbf,2)*0.01;

% Boxplot of the time-integrated number of leaders entering dangerous disks.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman');
boxchart(2*ones(size(nl_inside_Rk_mean2)), nl_inside_Rk_mean2);          % x = 1
boxchart(ones(size(nl_inside_nocbf_mean2)), nl_inside_nocbf_mean2);        % x = 2
xticks([1 2])
xticklabels({'NO CBF', 'CBF'})
ylabel('$\int_T N_L^{\mathrm{in}} \, dt$', 'Interpreter', 'latex')

% Boxplot of the time-integrated number of followers entering dangerous disks.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
set(gca, 'FontSize', 10, 'FontName', 'Times New Roman');
boxchart(2*ones(size(nf_inside_Rk_mean2)), nf_inside_Rk_mean2);          % x = 1
boxchart(ones(size(nf_inside_nocbf_mean2)), nf_inside_nocbf_mean2);        % x = 2
xticks([1 2])
xticklabels({'NO CBF', 'CBF'})
ylabel('$\int_T N_F^{\mathrm{in}} \, dt$', 'Interpreter', 'latex')


nf_inside_Rk_std = std(nf_inside_Rk, 0, 1) / 720;
nl_inside_Rk_std = std(nl_inside_Rk, 0, 1) / 100;
kl_leader_Rk_std = std(kl_leader_Rk, 0, 1);
kl_follower_Rk_std = std(kl_follower_Rk, 0, 1);
nf_inside_nocbf_std = std(nf_inside_nocbf, 0, 1) / 720;
nl_inside_nocbf_std = std(nl_inside_nocbf, 0, 1) / 100;
kl_leader_nocbf_std = std(kl_leader_nocbf, 0, 1);
kl_follower_nocbf_std = std(kl_follower_nocbf, 0, 1);

t = 0:0.01:100;

%% Follower
% Mean +/- one standard deviation of the follower dangerous-disk incursion rate.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
% Ensure proper scaling for LaTeX documents
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
% Configure tick labels with LaTeX interpreter and fontsize
set(gca, 'FontSize', 10, 'FontName', "Times New Roman");
h1 = fill([t fliplr(t)], ...
     [nf_inside_nocbf_mean + nf_inside_nocbf_std, fliplr(nf_inside_nocbf_mean - nf_inside_nocbf_std)], ...
     [0.0660 0.4430 0.7450], 'EdgeColor', [0.0660 0.4430 0.7450], 'FaceAlpha', 0.3);
plot(t, nf_inside_nocbf_mean, 'Color', '[0.0660 0.4430 0.7450]', 'LineWidth', 1.5);
h2 = fill([t fliplr(t)], ...
     [nf_inside_Rk_mean + nf_inside_Rk_std, fliplr(nf_inside_Rk_mean - nf_inside_Rk_std)], ...
     [0.8510 0.3255 0.0980], 'EdgeColor', [0.8510 0.3255 0.0980], 'FaceAlpha', 0.3);
plot(t, nf_inside_Rk_mean, 'Color', [0.8510 0.3255 0.0980], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\langle N_F^{\mathrm{in}}\rangle / N_F$', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');
lgd = legend([h1, h2], {'NO CBF', 'CBF'},  'Orientation','horizontal', ...
    'Location','northoutside');
lgd.FontSize = 8;
lgd.ItemTokenSize = [8 6];


% Mean +/- one standard deviation of the follower KL divergence.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
% Ensure proper scaling for LaTeX documents
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
% Configure tick labels with LaTeX interpreter and fontsize
set(gca, 'FontSize', 10, 'FontName', "Times New Roman");
h1 = fill([t fliplr(t)], ...
     [kl_follower_nocbf_mean + kl_follower_nocbf_std, fliplr(kl_follower_nocbf_mean - kl_follower_nocbf_std)], ...
     [0.0660 0.4430 0.7450], 'EdgeColor', [0.0660 0.4430 0.7450], 'FaceAlpha', 0.3);
plot(t, kl_follower_nocbf_mean, 'Color', '[0.0660 0.4430 0.7450]', 'LineWidth', 1.5);
h2 = fill([t fliplr(t)], ...
     [kl_follower_Rk_mean + kl_follower_Rk_std, fliplr(kl_follower_Rk_mean - kl_follower_Rk_std)], ...
     [0.8510 0.3255 0.0980], 'EdgeColor', [0.8510 0.3255 0.0980], 'FaceAlpha', 0.3);
plot(t, kl_follower_Rk_mean, 'Color', [0.8510 0.3255 0.0980], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$D_{KL}(\rho_F||\bar{\rho}_F)$', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');
lgd = legend([h1, h2], {'NO CBF', 'CBF'},  'Orientation','horizontal', ...
    'Location','northoutside');
lgd.FontSize = 8;
lgd.ItemTokenSize = [8 6];



%% Leader
% Mean +/- one standard deviation of the leader dangerous-disk incursion rate.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
% Ensure proper scaling for LaTeX documents
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
% Configure tick labels with LaTeX interpreter and fontsize
set(gca, 'FontSize', 10, 'FontName', "Times New Roman");
h1 = fill([t fliplr(t)], ...
     [nl_inside_nocbf_mean + nl_inside_nocbf_std, fliplr(nl_inside_nocbf_mean - nl_inside_nocbf_std)], ...
     [0.0660 0.4430 0.7450], 'EdgeColor', [0.0660 0.4430 0.7450], 'FaceAlpha', 0.3);
plot(t, nl_inside_nocbf_mean, 'Color', '[0.0660 0.4430 0.7450]', 'LineWidth', 1.5);
h2 = fill([t fliplr(t)], ...
     [nl_inside_Rk_mean + nl_inside_Rk_std, fliplr(nl_inside_Rk_mean - nl_inside_Rk_std)], ...
     [0.8510 0.3255 0.0980], 'EdgeColor', [0.8510 0.3255 0.0980], 'FaceAlpha', 0.3);
plot(t, nl_inside_Rk_mean, 'Color', [0.8510 0.3255 0.0980], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\langle N_L^{\mathrm{in}}\rangle / N_L$', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');
lgd = legend([h1, h2], {'NO CBF', 'CBF'},  'Orientation','horizontal', ...
    'Location','northoutside');
lgd.FontSize = 8;
lgd.ItemTokenSize = [8 6];


% Mean +/- one standard deviation of the leader KL divergence.
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
% Ensure proper scaling for LaTeX documents
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
% Configure tick labels with LaTeX interpreter and fontsize
set(gca, 'FontSize', 10, 'FontName', "Times New Roman");
h1 = fill([t fliplr(t)], ...
     [kl_leader_nocbf_mean + kl_leader_nocbf_std, fliplr(kl_leader_nocbf_mean - kl_leader_nocbf_std)], ...
     [0.0660 0.4430 0.7450], 'EdgeColor', [0.0660 0.4430 0.7450], 'FaceAlpha', 0.3);
plot(t, kl_leader_nocbf_mean, 'Color', '[0.0660 0.4430 0.7450]', 'LineWidth', 1.5);
h2 = fill([t fliplr(t)], ...
     [kl_leader_Rk_mean + kl_leader_Rk_std, fliplr(kl_leader_Rk_mean - kl_leader_Rk_std)], ...
     [0.8510 0.3255 0.0980], 'EdgeColor', [0.8510 0.3255 0.0980], 'FaceAlpha', 0.3);
plot(t, kl_leader_Rk_mean, 'Color', [0.8510 0.3255 0.0980], 'LineWidth', 1.5);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$D_{KL}(\rho_L||\bar{\rho}_L)$', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');
lgd = legend([h1, h2], {'NO CBF', 'CBF'},  'Orientation','horizontal', ...
    'Location','northoutside');
lgd.FontSize = 8;
lgd.ItemTokenSize = [8 6];

column_width = 0.75*7*0.49; % Column width in inches
aspect_ratio = 0.6;     % height / width
figure_height = column_width * aspect_ratio;

% Fraction of followers inside the target region.
chi_std = std(chi, 0, 1);
chi_mean = mean(chi,1);
chi_nocbf_std = std(chi_nocbf, 0, 1);
chi_nocbf_mean = mean(chi_nocbf,1);
figure('Units', 'inches', 'Position', [1, 1, column_width, figure_height]);
% Ensure proper scaling for LaTeX documents
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0, 0, column_width, figure_height]);
set(gcf, 'PaperSize', [column_width, figure_height]);
hold on;
% Configure tick labels with LaTeX interpreter and fontsize
set(gca, 'FontSize', 10, 'FontName', "Times New Roman");
h1 = fill([t fliplr(t)], ...
     [chi_nocbf_mean + chi_nocbf_std, fliplr(chi_nocbf_mean - chi_nocbf_std)], ...
     [0.0660 0.4430 0.7450], 'EdgeColor', [0.0660 0.4430 0.7450], 'FaceAlpha', 0.3);
plot(t, chi_nocbf_mean, 'Color', [0.0660 0.4430 0.7450], 'LineWidth', 1.2);
h2 = fill([t fliplr(t)], ...
     [chi_mean + chi_std, fliplr(chi_mean - chi_std)], ...
     [0.8510 0.3255 0.0980], 'EdgeColor', [0.8510 0.3255 0.0980], 'FaceAlpha', 0.3);
plot(t, chi_mean, 'Color', [0.8510 0.3255 0.0980], 'LineWidth', 1.2);
xlabel('$t$', 'Interpreter', 'latex');
ylabel('$\chi$', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter','latex');
lgd = legend([h1, h2], {'NO CBF', 'CBF'},  'Orientation','horizontal', ...
    'Location','southeast');
lgd.FontSize = 8;
% lgd = legend([h1, h2], {'NO CBF', 'CBF'}, 'Location', 'northwest');
lgd.ItemTokenSize = [8 6];

%% Export all figures to PDF

% Export every currently open figure in vector format for the manuscript.
fig_folder = fullfile(pwd,'figure');

if ~exist(fig_folder,'dir')
    mkdir(fig_folder);
end

figs = findall(groot,'Type','figure');

for k = 1:length(figs)

    fig = figs(k);

    filename = fullfile(fig_folder, sprintf('figure_%02d.pdf',k));

    exportgraphics(fig, filename, ...
        'ContentType','vector', ...
        'BackgroundColor','none');

end

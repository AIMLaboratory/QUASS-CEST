% plot barchart based on Delta Z output
%   20250617 Huabin Zhang, huabinz@connect.hku.hk
%       compare amide, guan, NOE
%       compare MTR, LDA, DMPLF
% clear
% close all

%% calculate contribution and sensitivity
list = ["output_MTR.mat","output_LDA.mat","output_DMPLF.mat",...
        "output_MTR_QUASS.mat","output_LDA_QUASS.mat","output_DMPLF_QUASS.mat"];
folder = ".\ran5000-B1_0.5_Ts_1\"; figName = "fig3_contri_sens_Ts1";
% folder = ".\ran5000-B1_0.5_Ts_2\"; figName = "figS1_contri_sens_Ts2";

tag = {'MTR_{asym}','LDA','DMPLF','Q-MTR_{asym}','Q-LDA','Q-DMPLF'};

contri_amide_all = [];
contri_guan_all = [];
contri_noe_all = [];
sens_amide_all = zeros(length(list),2); % slope and std
sens_guan_all = zeros(length(list),2);
sens_noe_all = zeros(length(list),2);
for idx = 1:length(list)
    load(folder+list(idx),'model_4var_amide','model_4var_guan','model_4var_noe');
    
    % contribution based on ANOVA
    anovaTable = anova(model_4var_amide);
    contri_amide_all = [contri_amide_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)]; % one column for one algorithm

    anovaTable = anova(model_4var_guan);
    contri_guan_all = [contri_guan_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)];

    anovaTable = anova(model_4var_noe);
    contri_noe_all = [contri_noe_all, 100 * anovaTable.SumSq(1:4) / sum(anovaTable.SumSq)];

    % sensitivity based on linear regression slope, unit [%/mM]
    DeltaZ_guan = model_4var_guan.Variables{:,end};
    DeltaZ_amide = model_4var_amide.Variables{:,end};
    DeltaZ_noe = model_4var_noe.Variables{:,end};

    linFit = fitlm(model_4var_guan.Variables{:,1},DeltaZ_guan);
    sens_guan_all(idx,1) = linFit.Coefficients.Estimate(2); % [%/mM]
    sens_guan_all(idx,2) = linFit.Coefficients.SE(2);

    linFit = fitlm(model_4var_amide.Variables{:,2},DeltaZ_amide);
    sens_amide_all(idx,1) = linFit.Coefficients.Estimate(2);
    sens_amide_all(idx,2) = linFit.Coefficients.SE(2);

    linFit = fitlm(model_4var_noe.Variables{:,4},DeltaZ_noe);
    sens_noe_all(idx,1) = linFit.Coefficients.Estimate(2);
    sens_noe_all(idx,2) = linFit.Coefficients.SE(2);
    
end

%% draw contribution and sensitivity
barthick = 0.6;
legendfontsize = 12;
colorset=[[085,059,148];[152,114,202];[246,231,237];[216,163,152]];
ylim_contri = [0,140];
ylim_sen = [0.001,0.04]; % [%/mM]

Fig2 = figure();set(gcf,'Position',[150 550 1200 450]);
tiledlayout(1,3,"TileSpacing","loose","Padding","loose")
nexttile

% 1-1 amide contribution
contri_amide_show = contri_amide_all([2,1,3,4],:); % resort display order
colorset_amide = colorset([2,1,3,4],:);

yyaxis left;
a = bar(1:length(list),contri_amide_show',barthick, 'stacked'); 
ylim(ylim_contri)
% title('Amide @ 3.5 ppm');
set(gca, 'xticklabel', tag,'FontSize',10,'FontWeight','bold');
ylabel('Variance Contribution [%]','FontSize',15);

for i = 1:length(a)
%     set(a(i),'FaceColor',colorset_amide(i,:)/255,'edgecolor','none');
    set(a(i),'FaceColor',colorset_amide(i,:)/255);
end

% 1-2 amide sensitivity
hold on;
yyaxis right;
plot(1:length(list), sens_amide_all(:,1), '-o', 'LineWidth', 2);
ylim(ylim_sen);

lgd = legend({'Amide', 'CEST@2ppm', 'MT', 'rNOE', 'Sensitivity'}, 'Location', 'northeast','Orientation','vertical','FontSize',legendfontsize); 

% 2-1 guan contribution
nexttile
yyaxis left;
a = bar(1:length(list),contri_guan_all',barthick, 'stacked'); 
ylim(ylim_contri)
% title('Guan @ 2.0 ppm');
set(gca, 'xticklabel', tag,'FontSize',10,'FontWeight','bold');

for i = 1:length(a)
    set(a(i),'FaceColor',colorset(i,:)/255);
end

% 2-2 guan sensitivity
hold on;
yyaxis right;
plot(1:length(list), sens_guan_all(:,1), '-o', 'LineWidth', 2);
ylim(ylim_sen);

% 3-1 noe contribution
nexttile
contri_noe_show = contri_noe_all([4,2,3,1],:); % resort display order
colorset_noe = colorset([4,2,3,1],:);

yyaxis left;
a = bar(1:length(list),contri_noe_show',barthick, 'stacked'); 
ylim(ylim_contri)
% title('NOE @ -3.5ppm');
set(gca, 'xticklabel', tag,'FontSize',10,'FontWeight','bold');

for i = 1:length(a)
    set(a(i),'FaceColor',colorset_noe(i,:)/255);
end

% 3-2 noe sensitivity
hold on;
yyaxis right;
plot(1:length(list), sens_noe_all(:,1), '-o', 'LineWidth', 2);
ylabel('Sensitivity [%/mM]','FontSize',15);
ylim(ylim_sen);


%% scatter plot
Fig1 = figure();set(gcf,'Position',[150 50 1200 450]);
tiledlayout('flow',"TileSpacing","compact","Padding","compact")
for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_amide');
    mdl = model_4var_amide;

    legend_on = 0;
    titletag = tag(i);
    xtag = 'f: Amide [mM]';
    ytag = '\Delta Z [%]';

    if i == 2 || i == 5 % LDA
        R2flag = 0; % average positive and negative respectively
    else
        R2flag = 1;
    end
    linearRegScatterPlot(mdl.Variables{:,2},mdl.Variables{:,end},xtag,ytag,titletag,legend_on,R2flag);
end

for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_guan');
    mdl = model_4var_guan;

    legend_on = 0;
    titletag = tag(i);
    xtag = 'f: CEST@2ppm [mM]';
    ytag = '\Delta Z [%]';

    if i == 2 || i == 5 % LDA
        R2flag = 0; % average positive and negative respectively
    else
        R2flag = 1;
    end
    linearRegScatterPlot(mdl.Variables{:,1},mdl.Variables{:,end},xtag,ytag,titletag,legend_on,R2flag);
end

for i = 1:length(list)
    nexttile
    load(folder+list(i),'model_4var_noe');
    mdl = model_4var_noe;

    legend_on = (i==length(list));
    titletag = tag(i);
    xtag = 'f: rNOE [mM]';
    ytag = '\Delta Z [%]';

    if i == 2 || i == 5 % LDA
        R2flag = 0; % average positive and negative respectively
    else
        R2flag = 1;
    end
    linearRegScatterPlot(mdl.Variables{:,4},mdl.Variables{:,end},xtag,ytag,titletag,legend_on,R2flag);
end

%% save
resol = '600';
OutputDir = ".\Out\";

if ~exist(DirTemp,'dir')
    mkdir(OutputDir);
end

fprintf("save to folder: "+join( split(OutputDir,'\'),'\\')+"\n");
exportgraphics(Fig1, OutputDir+figName+"_scatterplot.png", 'BackgroundColor', 'white', 'Resolution', resol);
exportgraphics(Fig2, OutputDir+figName+"_contribution.png", 'BackgroundColor', 'white', 'Resolution', resol);

%%
function p = linearRegScatterPlot(X,Y,xtag,ytag,titletag,legend_on,R2flag)
    % linear regression
    linFit = fitlm(X, Y);
    Y_fitted = linFit.Fitted;
   

    % R2 calculation
    warning off
    if ~isempty(find(Y<0, 1)) && R2flag == 0
        % for LDA, the reference signal would be influenced by MT effect,
        % which is prominent when B1 is high, resulting negative Lorentzian
        % difference. In such case, a seperate R2 calculation process is
        % employed.
        Yposidx = find(Y>=0);
        Ynegidx = find(Y<0);
    
        linFit_pos = fitlm(X(Yposidx),Y(Yposidx));
        linFit_neg = fitlm(X(Ynegidx),Y(Ynegidx));
    
        R2_pos = linFit_pos.Rsquared.Ordinary;
        R2_neg = linFit_neg.Rsquared.Ordinary;

        if isnan(R2_pos)
            R2_pos = 0;
        end
        if isnan(R2_neg)
            R2_neg = 0;
        end
        R2 = R2_pos * length(Yposidx)/length(Y) + R2_neg * length(Ynegidx)/length(Y);
    else
        linFit = fitlm(X, Y);
        R2 = linFit.Rsquared.Ordinary;
    end
    
    % sort data for plotting
    [X1_sorted, idx] = sort(X);
    Y_fitted_sorted = Y_fitted(idx);


    % scatter plot without CB
    scatter(X, Y, 2, 'b', 'filled'); hold on;
    p = plot(X1_sorted, Y_fitted_sorted, 'r-', 'LineWidth', 2);
    xlabel(xtag,'FontWeight','bold'); ylabel(ytag,'FontWeight','bold'); 
    % title(titletag);
    text(0.05, 0.93, ['R^2 = ' num2str(R2, '%.3f')], 'Units', 'normalized', 'Color', 'red', 'FontSize', 10,'FontWeight','bold');
    if legend_on
        lgd = legend('data', 'Fitted','Location','southeast');
        lgd.Position(1) = lgd.Position(1) + 0.01;
        lgd.Position(2) = lgd.Position(2) - 0.02;
    end
    hold off;
end




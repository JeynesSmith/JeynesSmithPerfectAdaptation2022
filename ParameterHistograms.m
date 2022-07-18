% The script is used to generate histograms for Figures 11 (KB1,KB2), and
% 12 (kA1, kA2, KA1, KA2) of the manuscript, and other parameter 
% histograms not used in the manuscript. This script loads the 
% datasets, determines which parameter sets are associated with RPA vs 
% non-RPA, and further splits RPA parameter sets into 3 bins based on the 
% magnitude of the RPA Range and the bounds defined by 'Bounds'. 

%% Load data and plot
% load in the dataset, define parameters of interest, a string for parameter names, x limits for histograms, and Bounds used for binning the RPA Ranges
load CCRSKms, SelectedParas = 5:8; LegString = {'K_{A1}','K_{A2}','K_{B1}','K_{B2}'}; XLimit = [-3,4];Bounds = [1,-2]; % michaelis constants
% load CCRSCat, SelectedParas = 1:4; LegString = {'k_{A1}','k_{A2}','k_{B1}','k_{B2}'}; XLimit = [-3,4];Bounds = [1,-3]; % catalytic constants
% load CCRSTot, SelectedParas = [9 11 10]; LegString = {'A_{tot}','E_{1tot}','B_{tot}'}; XLimit = [0,4];Bounds = [1,-2]; % total abundances

% define colours for plots
colours = [0 0.4470 0.7410;0.85 0.325 0.098;0.292 0.694 0.125;0.494 0.184 0.556]; linewidths = [2,1,2,1]; plotalphas = [0.8,0.4,0.8,0.4];

% Calculate the RPA Range
Range = RPASets(:,end) - RPASets(:,end-1);

% determine RPA and nonRPA psets
RPASets = RPASets(Range>0,SelectedParas); 
nonRPASets = nonRPASets(Range<=0,SelectedParas);
Range = log10(Range(Range>0));

% determine Bounds based on histogram
% figure
% histogram(Range,'binwidth',1)
% keyboard

% plot RPA vs nonRPA
f = figure;
f.Position = [38.6,211.4,741.6,541.6];
t = tiledlayout(2,2,'TileSpacing','Compact');
% RPA 
for i = 1:length(SelectedParas)
    if i == 1
        nexttile(1)
        hold on
    elseif i == 3
        nexttile(3)
        hold on
    end
    histogram(log10(RPASets(:,i)),'binwidth',1,'FaceColor',colours(i,:),'EdgeColor',colours(i,:),'FaceAlpha',plotalphas(i),'LineWidth',linewidths(i))
    if i == 2
        legend(LegString(1:2),'fontsize',12),ylim([0,length(RPASets)/5*2])
        xlim(XLimit)
        title('RPA Behaviour','fontsize',12,'fontweight','normal')
    elseif i == length(SelectedParas)
        legend(LegString(3:end),'fontsize',12),ylim([0,length(RPASets)/5*2])
        xlim(XLimit)
    end
end
% nonRPA 
for i = 1:length(SelectedParas)
    if i == 1
        nexttile(2)
        hold on
    elseif i == 3
        nexttile(4)
        hold on
    end
    histogram(log10(nonRPASets(:,i)),'binwidth',1,'FaceColor',colours(i,:),'EdgeColor',colours(i,:),'FaceAlpha',plotalphas(i),'LineWidth',linewidths(i))
    if i == 2
        title('No RPA Behaviour','fontsize',11,'fontweight','normal')
        legend(LegString(1:2),'fontsize',12),ylim([0,length(nonRPASets)/5*2]),xlim(XLimit)
    elseif i == length(SelectedParas)
        legend(LegString(3:end),'fontsize',12),ylim([0,length(nonRPASets)/5*2]),xlim(XLimit)
    end
end
ylabel(t,'Occurence','fontsize',15)
xlabel(t,'log(Constants)','fontsize',15)


% plot RPA ranges
f = figure;
f.Position = [38.6,211.4,741.6,541.6];
t = tiledlayout(2,3,'TileSpacing','Compact');
Bounds = [Inf Bounds -Inf];
for k = 1:2
    for j = 1:3
        nexttile((k-1)*3+j)
        hold on
        for i = (k-1)*2+(1:2)
            if length(SelectedParas)==3 && i==4, continue, end
            TempVec = RPASets(Range<Bounds(j) & Range>=Bounds(j+1),i);
            histogram(log10(TempVec),'binwidth',1,'FaceColor',colours(i,:),'EdgeColor',colours(i,:),'FaceAlpha',plotalphas(i),'LineWidth',linewidths(i))
        end
        if k==1,legend(LegString(1:2),'fontsize',12),end
        if k==2,legend(LegString(3:end),'fontsize',12),end
        xlim(XLimit),ylim([0,max(length(TempVec),1)])
        if j==1&&k==1,title(['RPA Range > 10^{' num2str(Bounds(j+1)) '}'],'fontsize',11,'fontweight','normal'),end
        if j==2&&k==1,title(['10^{' num2str(Bounds(j)) '} > RPA Range > 10^{' num2str(Bounds(j+1)) '}'],'fontsize',11,'fontweight','normal'),end
        if j==3&&k==1,title(['RPA Range < 10^{' num2str(Bounds(j)) '}'],'fontsize',11,'fontweight','normal'),end
    end
end
xlabel(t,'log_{10}(Parameters)','fontsize',15)
ylabel(t,'Occurence','fontsize',15)

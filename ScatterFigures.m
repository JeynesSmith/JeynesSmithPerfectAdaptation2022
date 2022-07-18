% The script is used to generate scatter plots for Figures 7 (Ma RPA), 
% 8 (CC RPA), and 10 (MA vs CC RPA Range) of the manuscript. This script 
% loads the datasets used to examine the roles of total substrate and
% enzyme abundances on RPA presence and the RPA range. The magnitude of the
% range is represented by the colour, defined in the following section.

%% colors
customcolor1 = [240,249,232
204,235,197
168,221,181
123,204,196
78,179,211
43,140,190
8,88,158]./255;
customcolor2 = [254,240,217
253,212,158
253,187,132
252,141,89
227,74,51
179,0,0]./255;
customcolor3 = [255,255,178
254,204,92
253,141,60
240,59,32
189,0,38]./255;

%% figures
% complex complete: low sensitivity 
load CCRSTot RPASets
CCRange = RPASets(:,end) - RPASets(:,end-1); % find the range
CCRPASets = RPASets(CCRange>0,:); % find appropriate datasets with positive range
Range = CCRange(CCRange>0);
figure % Atot vs Etot
scatter(log10(CCRPASets(:,9)),log10(CCRPASets(:,11)),20,min([log10(Range), 5.*ones(size(Range))],[],2),'filled')
colorbar
caxis([-2,5]),colormap(customcolor1)
xlabel('log_{10}(A_{tot})'),ylabel('log_{10}(E_{tot})')
figure % Btot vs Etot
scatter(log10(CCRPASets(:,10)),log10(CCRPASets(:,11)),20,min([log10(Range), 5.*ones(size(Range))],[],2),'filled')
colorbar
caxis([-2,5]),colormap(customcolor1)
xlabel('log_{10}(B_{tot})'),ylabel('log_{10}(E_{tot})')

% Ma et al.: low sensitivity 
load MaRSTot RPASets
MaRange = RPASets(:,end) - RPASets(:,end-1); % find the range
MaRPASets = RPASets(MaRange>0,:); % find appropriate datasets with positive range
Range = MaRange(MaRange>0);
figure % Atot vs Etot
scatter(log10(MaRPASets(:,9)),log10(MaRPASets(:,11)),20,min([log10(Range), 5.*ones(size(Range))],[],2),'filled')
colorbar
caxis([-1,5]),colormap(customcolor2)
xlabel('log_{10}(A_{tot})'),ylabel('log_{10}(E_{tot})')
figure % Btot vs Etot
scatter(log10(MaRPASets(:,10)),log10(MaRPASets(:,11)),20,min([log10(Range), 5.*ones(size(Range))],[],2),'filled')
colorbar
caxis([-1,5]),colormap(customcolor2)
xlabel('log_{10}(B_{tot})'),ylabel('log_{10}(E_{tot})')

% range comparison for low sensitivity 
load CCRSTot RPASets
SharedParas = RPASets(MaRange>0 & CCRange>0,:); % find parameter sets where the range is positive in Ma and CC models
SharedRange = CCRange(MaRange>0 & CCRange>0)./MaRange(MaRange>0 & CCRange>0); % find positive ranges where the range is positive in Ma and CC models
SharedRange = max([SharedRange,1./SharedRange],[],2); % set the ranges to be always greater than 1 (i.e. large range / small range)
figure % Etot/Atot vs Etot/Btot
scatter(log10(SharedParas(:,11)./SharedParas(:,9)),log10(SharedParas(:,11)./SharedParas(:,10)),20,min([log10(SharedRange), 5.*ones(size(SharedRange))],[],2),'filled')
colorbar
caxis([0,5]),colormap(customcolor3)
xlabel('log_{10}(E_{tot}/A_{tot})'),ylabel('log_{10}(E_{tot}/B_{tot})')
ylim([-4,0])

% complex complete: high sensitivity 
load CCRSTot_SmallKm RPASets
CCRange = RPASets(:,end) - RPASets(:,end-1); % find the range
CCRPASets = RPASets(CCRange>0,:); % find appropriate datasets with positive range
CCRange = CCRange(CCRange>0);
figure % Atot vs Etot
scatter(log10(CCRPASets(:,9)),log10(CCRPASets(:,11)),20,min([log10(CCRange), 5.*ones(size(CCRange))],[],2),'filled')
colorbar
caxis([-2,5]),colormap(customcolor1)
xlabel('log_{10}(A_{tot})'),ylabel('log_{10}(E_{tot})')
figure % Btot vs Etot
scatter(log10(CCRPASets(:,10)),log10(CCRPASets(:,11)),20,min([log10(CCRange), 5.*ones(size(CCRange))],[],2),'filled')
colorbar
caxis([-2,5]),colormap(customcolor1)
xlabel('log_{10}(B_{tot})'),ylabel('log_{10}(E_{tot})')

% Ma et al.: high sensitivity 
load MaRSTot_SmallKm RPASets
MaRange = RPASets(:,end) - RPASets(:,end-1); % find the range
MaRPASets = RPASets(MaRange>0,:); % find appropriate datasets with positive range
MaRange = MaRange(MaRange>0);
figure % Atot vs Etot
scatter(log10(MaRPASets(:,9)),log10(MaRPASets(:,11)),20,min([log10(MaRange), 5.*ones(size(MaRange))],[],2),'filled')
colorbar
caxis([-1,5]),colormap(customcolor2)
xlabel('log_{10}(A_{tot})'),ylabel('log_{10}(E_{tot})')
figure
scatter(log10(MaRPASets(:,10)),log10(MaRPASets(:,11)),20,min([log10(MaRange), 5.*ones(size(MaRange))],[],2),'filled')
colorbar % Btot vs Etot
caxis([-1,5]),colormap(customcolor2)
xlabel('log_{10}(B_{tot})'),ylabel('log_{10}(E_{tot})')

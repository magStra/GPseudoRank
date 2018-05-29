%look at precocious cells
%their indices are 91 and 92
load('ShalekUncert.mat');%we use 500,000 samples from one MCMC chain
%produced e.g. by ShalekFewerCaptureTimes.m
xx(91,:);
xx(92,:);
%plot the posterior distributions of the precicous cells identified in 
%Shalek et al. Nature. 2014
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 8)
set(0,'defaultfigurecolor',[1 1 1]);
figure();
bar(1:307,posFreq(91,:,1),'FaceAlpha',0.6,'FaceColor','g')
hold on;
bar(1:307,posFreq(92,:,1),'FaceAlpha',0.6,'FaceColor','b')
hold on;
a=ylim;
h1 = line([125,125],[0,a(2)]);
h1.LineWidth = 2;
h1.Color = 'm';
xlabel('pseudo-position','FontSize',10);
ylabel('frequency','FontSize',10);
leg = legend('1h: S51', '1h: S52','cap. time 2h','Location','NorthOutSide',...
    'Orientation','horizontal');
legend boxoff
xlim([80,200]);
leg.FontSize = 10;
box off;
set(gcf, 'PaperUnits', 'centimeters');
 x_width=8.8 ;y_width=5.5;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 set(gcf,'PaperSize',[x_width y_width]);
print('PreCells','-dpdf','-r300');

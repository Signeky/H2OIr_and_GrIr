%Load GrIr(111) 125 K data.. 
clear all; close all; 

%% Load and fit data 

x = linspace(2.8,800,800);

addpath Dynamics\Ir_and_GrIr
fileListGrIr125 = {'dy019379.mat','dy019380.mat','dy019381.mat','dy019382.mat','dy019383.mat','dy019384.mat','dy019385.mat','dy019386.mat','dy019387.mat','dy019388.mat','dy019389.mat','dy019390.mat','dy019391.mat','dy019392.mat','dy019393.mat','dy019394.mat','dy019395.mat','dy019396.mat','dy019397.mat','dy019398.mat'};  

%Sort dKs from low to high 
for i = 1:length(fileListGrIr125)
    load(fileListGrIr125{i})
    dKGrIr125 = abs(meas.dK);
    dKsGrIr125(i) = dKGrIr125; 
  
end
[dKsGrIr125_sorted, dKsGrIr125_order] = sort(dKsGrIr125);
newfileListGrIr125 = fileListGrIr125(dKsGrIr125_order); %dK ordered file list 


for i=1:length(newfileListGrIr125)
    
    % load in file from the list above
    load(newfileListGrIr125{i})
    dKGrIr125 = meas.dK; 
    dK(i) = meas.dK; 


    TGrIr125 = str2num(meas.endStatus.tSample);
    tseGrIr125 = meas.setime;
    PmagGrIr125 = meas.mean.Pmag;
    PimGrIr125 = meas.mean.Pimag;
    
    
    %We cutoff at n=53 in the fit to exclude datapoints arising from quasi-elastic scattering from phonons
    cutoffn = 53;
    [xData, yData] = prepareCurveData(tseGrIr125(cutoffn:end),PmagGrIr125(cutoffn:end));
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    %%Force the fit to decay to find upper limit 
    opts.Lower = [0 0 PmagGrIr125(end)/2]; 
    opts.StartPoint = [0.5 0.01 PmagGrIr125(end)/2];    
    opts.Upper = [1 0.5 PmagGrIr125(end)/2];
   
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.TolFun = 1e-08; 
    
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Rsquare_GrIr125(i)=gof.rsquare;
    
    ci = confint(fitresult,0.68);
    
    %Store fitresult in arrays
    alpha_GrIr125(i)=fitresult.b;
    dalpha_GrIr125(i)=abs(ci(1,2)-ci(2,2))/2;
    
    offsets_GrIr125(i)=fitresult.c;
    doffsets_GrIr125(i)=abs(ci(1,3)-ci(2,3))/2;
   
    amplitude_GrIr125(i)=fitresult.a;
    damplitude_GrIr125(i)=abs(ci(1,1)-ci(2,1))/2;
     
    %Plot ISF for each delta K 
    figure(2)
    subplot(4,5,i)
    semilogx(tseGrIr125, PmagGrIr125,'s','color','#0072bd','MarkerSize',5);
    title(['\DeltaK = ' num2str(dKGrIr125,2)]);
    xlim([0,800]); ylim([0.98*min(yData),max(yData)*1.02]);
    hold on;
    semilogx(x, fitresult.a*exp(-fitresult.b*x)+fitresult.c,'color','#045388')
    
    legend('off');
    grid on
    hold off
   
end

    
%% Figure summarising fit parameters
figure2=figure('units','centimeters','position',[5,3,30,15],'color','white','DefaultLineLineWidth',1.5);
set(figure2,'DefaultLineLineWidth',2); axes2 = axes('Parent',figure2,'LineWidth',1.5,'FontSize',16);
box(axes2,'on'); hold(axes2,'all');

subplot(2,2,1)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsGrIr125_sorted,alpha_GrIr125,'bo','color','#045388','MarkerSize',8);
errorbar(dKsGrIr125_sorted,alpha_GrIr125,dalpha_GrIr125,'color',[0.4 0.4 0.4],'LineStyle','None');
ylim([0 0.0005]);
grid on; 
xlabel('\DeltaK (Å^{-1})');  ylabel('\alpha (ps^{-1})');

subplot(2,2,2)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsGrIr125_sorted,amplitude_GrIr125,'bo','color','#045388','MarkerSize',8);
errorbar(dKsGrIr125_sorted,amplitude_GrIr125,damplitude_GrIr125,'color',[0.4 0.4 0.4],'LineStyle','None');
xlim([0 3.2]); 

xlabel('\DeltaK (Å^{-1})');  ylabel('Amplitude');
grid on;

subplot(2,2,3)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsGrIr125_sorted,offsets_GrIr125,'bo','color','#045388','MarkerSize',8);
errorbar(dKsGrIr125_sorted,offsets_GrIr125,doffsets_GrIr125,'color',[0.4 0.4 0.4],'LineStyle','None');
xlim([0 3.2]); 
xlabel('\DeltaK (Å^{-1})');  ylabel('offset');
grid on;

subplot(2,2,4)

set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsGrIr125_sorted,Rsquare_GrIr125,'bo','color','#045388','MarkerSize',8);
xlim([0 3.2]); ylim([0 1]);
xlabel('\DeltaK (Å^{-1})');  ylabel('R^2');
grid on;



clear all; close all; 

x = linspace(2.8,800,800);


addpath Dynamics\Ir_and_GrIr
fileListIr135 = {'dy019490.mat','dy019491.mat','dy019492.mat','dy019493.mat','dy019494.mat','dy019495.mat','dy019496.mat','dy019497.mat','dy019498.mat','dy019499.mat','dy019500.mat','dy019501.mat','dy019502.mat','dy019503'};  

%Sort dKs from low to high 
for i = 1:length(fileListIr135)
    load(fileListIr135{i})
    dKIr135 = abs(meas.dK);
    dKsIr135(i) = dKIr135; 
  
end
[dKsIr135_sorted, dKsIr135_order] = sort(dKsIr135);
newfileListIr135 = fileListIr135(dKsIr135_order); %dK ordered file list 


for i=1:length(newfileListIr135)
    
    % load in file from the list above
    load(newfileListIr135{i})
    dKIr135 = meas.dK; 
    TIr135 = str2num(meas.endStatus.tSample);
    tseIr135 = meas.setime;
    PmagIr135 = meas.mean.Pmag;
    PimIr135 = meas.mean.Pimag;
    
    
    %We cutoff at n=2 in the fit to exclude datapoints arising from quasi-elastic scattering from phonons
    [xData, yData] = prepareCurveData(tseIr135(2:end),PmagIr135(2:end));
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.StartPoint = [0.4 0.01 min(PmagIr135)];
    opts.Upper = [1 0.5 0.6];
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.TolFun = 1e-08; 
    
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Rsquare_Ir135(i)=gof.rsquare;
    
    ci = confint(fitresult,0.68);
    
    %Store fitresult in arrays
    alpha_Ir135(i)=fitresult.b;
    dalpha_Ir135(i)=abs(ci(1,2)-ci(2,2))/2;
    
    offsets_Ir135(i)=fitresult.c;
    doffsets_Ir135(i)=abs(ci(1,3)-ci(2,3))/2;
   
    amplitude_Ir135(i)=fitresult.a;
    damplitude_Ir135(i)=abs(ci(1,1)-ci(2,1))/2;
     
    %Plot ISF for each delta K 
    figure(2)
    subplot(4,4,i)
    semilogx(tseIr135, PmagIr135,'s','color','#35845a','MarkerSize',5);
    title(['\DeltaK = ' num2str(dKIr135,2)]);
    xlim([0,705]); ylim([0.98*min(yData),max(yData)*1.02]);
    hold on;
    semilogx(x, fitresult.a*exp(-fitresult.b*x)+fitresult.c,'color','#286444')
    
    legend('off');
    grid on
    hold off
   
end

    
%Figure summarising fit parameters
figure2=figure('units','centimeters','position',[5,3,30,15],'color','white','DefaultLineLineWidth',1.5);
set(figure2,'DefaultLineLineWidth',2); axes2 = axes('Parent',figure2,'LineWidth',1.5,'FontSize',16);
box(axes2,'on'); hold(axes2,'all');

subplot(2,2,1)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsIr135_sorted,alpha_Ir135,'bo','color','#286444','MarkerSize',8);
errorbar(dKsIr135_sorted,alpha_Ir135,dalpha_Ir135,'color',[0.4 0.4 0.4],'LineStyle','None');
ylim([-0.001 0.017]);
xlabel('\DeltaK (Å^{-1})');  ylabel('\alpha (ps^{-1})');
grid on;
hold on

subplot(2,2,2)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsIr135_sorted,amplitude_Ir135,'bo','color','#286444','MarkerSize',8);
errorbar(dKsIr135_sorted,amplitude_Ir135,damplitude_Ir135,'color',[0.4 0.4 0.4],'LineStyle','None');
xlim([0 3.2]); 

xlabel('\DeltaK (Å^{-1})');  ylabel('Amplitude');
grid on;

subplot(2,2,3)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsIr135_sorted,offsets_Ir135,'bo','color','#286444','MarkerSize',8);
errorbar(dKsIr135_sorted,offsets_Ir135,doffsets_Ir135,'color',[0.4 0.4 0.4],'LineStyle','None');
xlim([0 3.2]); 
xlabel('\DeltaK (Å^{-1})');  ylabel('offset');
grid on;

subplot(2,2,4)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsIr135_sorted,Rsquare_Ir135,'bo','color','#286444','MarkerSize',8);
xlim([0 3.2]); ylim([0 1]);
xlabel('\DeltaK (Å^{-1})');  ylabel('R^2');
grid on;


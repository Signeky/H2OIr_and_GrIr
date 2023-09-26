
clear all; close all; 

%% Load and fit data 

x = linspace(2.8,800,800);


addpath Dynamics\Ir_and_GrIr
fileListIr125 = {'dy019476.mat','dy019477.mat','dy019478.mat','dy019479.mat','dy019480.mat','dy019481.mat','dy019482.mat','dy019483.mat','dy019484.mat','dy019485.mat','dy019486.mat','dy019487.mat','dy019488.mat','dy019489'};  

%Sort dKs from low to high 
for i = 1:length(fileListIr125)
    load(fileListIr125{i})
    dKIr125 = abs(meas.dK);
    dKsIr125(i) = dKIr125; 
  
end
[dKsIr125_sorted, dKsIr125_order] = sort(dKsIr125);
newfileListIr125 = fileListIr125(dKsIr125_order); %dK ordered file list 


for i=1:length(newfileListIr125)
    
    % load in file from the list above
    load(newfileListIr125{i})
    dKIr125 = meas.dK; 
    dK(i) = meas.dK; 


    TIr125 = str2num(meas.endStatus.tSample);
    tseIr125 = meas.setime;
    PmagIr125 = meas.mean.Pmag;
    PimIr125 = meas.mean.Pimag;
    
    %We cutoff at n=2 in the fit to exclude datapoints arising from quasi-elastic scattering from phonons  
    [xData, yData] = prepareCurveData(tseIr125(2:end),PmagIr125(2:end));
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.StartPoint = [0.4 0.01 min(PmagIr125(end))];
    opts.Upper = [1 0.5 0.6];
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.TolFun = 1e-08; 
    

    
    
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Rsquare_Ir125(i)=gof.rsquare;
    
    ci = confint(fitresult,0.68);
    
    %Store fitresult in arrays
    alpha_Ir125(i)=fitresult.b;
    dalpha_Ir125(i)=abs(ci(1,2)-ci(2,2))/2;
    
    offsets_Ir125(i)=fitresult.c;
    doffsets_Ir125(i)=abs(ci(1,3)-ci(2,3))/2;
   
    amplitude_Ir125(i)=fitresult.a;
    damplitude_Ir125(i)=abs(ci(1,1)-ci(2,1))/2;
     
    %Plot ISF for each delta K 
    figure(2)
    subplot(4,4,i)
    semilogx(tseIr125, PmagIr125,'s','color','#35845a','MarkerSize',5);
    title(['\DeltaK = ' num2str(dKIr125,2)]);
    xlim([0,705]); ylim([0.98*min(yData),max(yData)*1.02]);
    hold on;
    semilogx(x, fitresult.a*exp(-fitresult.b*x)+fitresult.c,'color','#286444')
    
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
plot(dKsIr125_sorted,alpha_Ir125,'bo','color','#286444','MarkerSize',8);
errorbar(dKsIr125_sorted,alpha_Ir125,dalpha_Ir125,'color',[0.4 0.4 0.4],'LineStyle','None');
ylim([-0.001 0.017]);
xlim([0 3.2]);
grid on;
hold on


%Fit alphas to the CE model 
[deltaks_Ir125, alphas_Ir125] = prepareCurveData(dKsIr125_sorted,alpha_Ir125);

ft = fittype( ['2/(t)*( (sin(x*2.72./2*cosd(30)).^2)+(sin(x*2.72./2*cosd(90)).^2)+(sin(x*2.72./2*cosd(150)).^2)+(sin(x*2.72./2*cosd(210)).^2)+' ...
    '(sin(x*2.72./2*cosd(270)).^2)+(sin(x*2.72./2*cosd(330)).^2) )'], 'independent', 'x', 'dependent', 'y', 'coefficients',{'t'} );


opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [65];
opts.StartPoint = [1300];
opts.Upper = [2000];
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.TolFun = 1e-08;
opts.Weights = 1./dalpha_Ir125;
[fitresult_sin, gof_sin] = fit( deltaks_Ir125, alphas_Ir125, ft, opts )

ci_sin = confint(fitresult_sin,0.68);
x_dKs = linspace(0,3.2);


fsin =(2/(fitresult_sin.t)*((sin(x_dKs*2.72./2*cosd(30)).^2)+(sin(x_dKs*2.72./2*cosd(90)).^2)+(sin(x_dKs*2.72./2*cosd(150)).^2)+ ...
    (sin(x_dKs*2.72./2*cosd(210)).^2)+(sin(x_dKs*2.72./2*cosd(270)).^2)+(sin(x_dKs*2.72./2*cosd(330)).^2)));


plot(x_dKs,fsin,'color',["#daa711"])
set(gca,'LineWidth',1.7,'FontSize',12,'Box', 'on');
xlabel('$\Delta K$(\AA$^{-1}$)');  ylabel('$\alpha$ (ps$^{-1}$)');

%Calculation of the diffusion values
tau = fitresult_sin.t*1e-12; %in s

dtau = abs(ci_sin(1,1)-ci_sin(2,1))/2*1e-12;

l1 = 2.72*10^(-10); %distance to nearest neighbour in m
dl1 = 0.03*10^(-10); %take from diffraction scan
diffusion_cst = (1/(4*tau))*l1^2;
ddiffusion_cst = sqrt(((l1^2/(4*tau^2))*dtau)^2+(((-l1/(2*tau))^2*(dl1)^2))); %through error proporgation

fprintf(['resident time tau = (%2.1f ' char(177) ' %2.1f) ps'],tau*10^(12), dtau*10^(12))

fprintf(['\njump length = (%2.1f ' char(177) ' %2.1f) Å'],l1*10^(10), dl1*10^(10))

fprintf(['\nDiffusion constant = (%2.0f ' char(177) ' %2.0f) pm^2/s\n'],diffusion_cst*10^(12), ddiffusion_cst*10^(12))



subplot(2,2,2)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsIr125_sorted,amplitude_Ir125,'bo','color','#286444','MarkerSize',8);
errorbar(dKsIr125_sorted,amplitude_Ir125,damplitude_Ir125,'color',[0.4 0.4 0.4],'LineStyle','None');
xlim([0 3.2]); 

xlabel('\DeltaK (Å^{-1})');  ylabel('Amplitude');
grid on;

subplot(2,2,3)
set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsIr125_sorted,offsets_Ir125,'bo','color','#286444','MarkerSize',8);
errorbar(dKsIr125_sorted,offsets_Ir125,doffsets_Ir125,'color',[0.4 0.4 0.4],'LineStyle','None');
xlim([0 3.2]); 
xlabel('\DeltaK (Å^{-1})');  ylabel('offset');
grid on;

subplot(2,2,4)

set(gca,'LineWidth',2,'FontSize',14,'Layer','top','Box','on'); hold(gca,'all');
plot(dKsIr125_sorted,Rsquare_Ir125,'bo','color','#286444','MarkerSize',8);
xlim([0 3.2]); ylim([0 1]);
xlabel('\DeltaK (Å^{-1})');  ylabel('R^2');
grid on;


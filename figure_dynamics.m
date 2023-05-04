
clear all; close all; 


%%Load data Ir(111) 125 K 
%addpath U:\Documents\MATLAB\OpenAccess\Dynamics
addpath Dynamics
fileListIr125 = {'dy019476.mat','dy019477.mat','dy019478.mat','dy019479.mat','dy019480.mat','dy019481.mat','dy019482.mat','dy019483.mat','dy019484.mat','dy019485.mat','dy019486.mat','dy019487.mat','dy019488.mat','dy019489'};  


for i=1:length(fileListIr125)
    
    % load in file from the list above
    load(fileListIr125{i})
    dKIr125(i) = meas.dK; 
     
    TIr125 = str2num(meas.endStatus.tSample);
    tseIr125 = meas.setime;
    PmagIr125 = meas.mean.Pmag;
    PimIr125 = meas.mean.Pimag;
 
    %Fit data 
    [xData, yData] = prepareCurveData(tseIr125(2:end),PmagIr125(2:end));
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
%     opts.Lower = [0 0 0];
%     opts.StartPoint = [0.4 0.01 0];
%     opts.Upper = [1 0.4 PmagIr125(end)];
%     
    opts.Lower = [0 0 0];
    opts.StartPoint = [0.4 0.01 min(PmagIr125)];
    opts.Upper = [1 0.5 0.6];
    
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.TolFun = 1e-08; 
    [fitresultIr125, gof] = fit( xData, yData, ft, opts );
    Rsquare_Ir125(i)=gof.rsquare;
    
    ci = confint(fitresultIr125,0.68);
    
    %Store fitresult in arrays
    alpha_Ir125(i)=fitresultIr125.b;
    dalpha_Ir125(i)=abs(ci(1,2)-ci(2,2))/2;
    
    %Store data and fit paramters at dK = 1 
    if (dKIr125(i)<1.1) && (dKIr125(i)>0.9)
        tseIr125_dK1 = tseIr125;
        PmagIr125_dK1 = PmagIr125;
        dK1_Ir125 = dKIr125; 
        a_Ir125 = fitresultIr125.a; 
        b_Ir125 = fitresultIr125.b; 
        c_Ir125 = fitresultIr125.c; 
    end

   
end



%% Load data GrIr(111) 125 K 

fileListGrIr125 = {'dy019379.mat','dy019380.mat','dy019381.mat','dy019382.mat','dy019383.mat','dy019384.mat','dy019385.mat','dy019386.mat','dy019387.mat','dy019388.mat','dy019389.mat','dy019390.mat','dy019391.mat','dy019392.mat','dy019393.mat','dy019394.mat','dy019395.mat','dy019396.mat','dy019397.mat','dy019398.mat'};  

for i=1:length(fileListGrIr125)
    
    % load in file from the list above
    load(fileListGrIr125{i})
    dKGrIr125(i) = meas.dK; 
    
    TGrIr125 = str2num(meas.endStatus.tSample);
    tseGrIr125 = meas.setime;
    PmagGrIr125 = meas.mean.Pmag;
    PimGrIr125 = meas.mean.Pimag;
    
    %Fit data 
    cutoffn = 53; %Cut off phonon data
    [xData, yData] = prepareCurveData(tseGrIr125(cutoffn:end),PmagGrIr125(cutoffn:end));
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [0 0 PmagGrIr125(end)/2];
    opts.StartPoint = [0.5 0.01 PmagGrIr125(end)/2];    
    opts.Upper = [1 0.5 PmagGrIr125(end)/2];
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.TolFun = 1e-08; 
    [fitresultGrIr, gof] = fit( xData, yData, ft, opts );
    Rsquare_GrIr125(i)=gof.rsquare;
    
    ci = confint(fitresultGrIr,0.68);
    

    %Store fitresult in arrays
    alpha_GrIr125(i)=fitresultGrIr.b;
    dalpha_GrIr125(i)=abs(ci(1,2)-ci(2,2))/2;
    
    %Store data and fitresult for dK = 1
    if (dKGrIr125(i)<1.1) && (dKGrIr125(i)>0.9)
        tseGrIr125_dK1 = tseGrIr125;
        PmagGrIr125_dK1 = PmagGrIr125;
        dK1_GrIr125 = dKGrIr125; 
        a_GrIr125 = fitresultGrIr.a; 
        b_GrIr125 = fitresultGrIr.b; 
        c_GrIr125 = fitresultGrIr.c; 
    end
   
end



addpath Dynamics\GrNi

%%Load data GrNi(111) 
fileListGrNi125 = [10607:10620,10622:10626,10627:10632,10634:10638,10640,10641,10643,10644,10621];

for kk=1:length(fileListGrNi125)
    
    %=================== Load .mat file and extract information ===================
    filename=char(['dy0',num2str(fileListGrNi125(kk),'%05d')]);
    disp(['Loading #' num2str(fileListGrNi125(kk))]);
    load([filename,'.mat']);
    
    dK_GrNi125(kk) = abs(data.dK);
    

    cutoff(kk)=data.cutoff;
    xdata = data.xdata;
    ydata = data.ydata;
    
    normfact(kk)=ydata(1);
    ydata(end-4:end)=mean(ydata(end-4:end));
    
    if dK_GrNi125(kk)<=0.06;
       cutoffn=15;
    elseif (dK_GrNi125(kk) > 0.06) && (dK_GrNi125(kk) < 0.15);
       cutoffn=20;
    elseif (dK_GrNi125(kk) > 0.7) && (dK_GrNi125(kk) < 0.8);
       cutoffn=20;
    elseif (dK_GrNi125(kk) > 0.8) && (dK_GrNi125(kk) < 1.2);
       cutoffn=12;
    elseif (dK_GrNi125(kk) > 1.2) && (dK_GrNi125(kk) < 2.7);
       cutoffn=10;
    else
       cutoffn=12;
    end

    cutoffvec(kk)=cutoffn;
   
    
    x=xdata(xdata>cutoffn); y=ydata(xdata>cutoffn);
    [xData, yData] = prepareCurveData(x,y);

    % Set up fittype and options.
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.StartPoint = [0.5 0.05 min(y)];
    opts.Upper = [1 0.5 0.6];
    opts.MaxFunEvals = 1000;
    opts.MaxIter = 1000;
    opts.TolFun = 1e-08;
   
    [fitresult, gof] = fit( xData, yData, ft, opts );
    Rsquare(kk)=gof.rsquare;
 
    %Store data at dK = 1 
    if (dK_GrNi125(kk)<1) && (dK_GrNi125(kk)>0.98)
        index_GrNi125 = kk;
        tseGrNi125_dK1 = xdata;
        PmagGrNi125_dK1 = ydata;
        dK1_GrNi125 = dK_GrNi125(kk);
    
        
    end
    
    ci = confint(fitresult,0.68);
    
    %Store fit parametrs
    Amp_GrNi125(kk)=fitresult.a;
    dAmp_GrNi125(kk)=abs(ci(1,1)-ci(2,1))/2;
    alpha_GrNi125(kk)=fitresult.b;
    dalpha_GrNi125(kk)=abs(ci(1,2)-ci(2,2))/2;
    offset_GrNi125(kk)=fitresult.c;
    doffset_GrNi125(kk)=abs(ci(1,3)-ci(2,3))/2;
    
    
end

%Create x array when plotting fits 
x = linspace(8,800);
x_dKs = linspace(0,3.2);

figure=figure('units','centimeters','position',[0,0,32,16],'color','white','DefaultTextInterpreter','LaTex','DefaultLineLineWidth',1.5);
set(figure,'DefaultLineLineWidth',2); axes2 = axes('Parent',figure,'LineWidth',1.5,'FontSize',14);

subplot(2,3,1)
IrCurve = semilogx(tseIr125_dK1,PmagIr125_dK1,'s','color',[0.9290 0.6940 0.1250],'displayname',['Ir(111), \DeltaK = ' num2str(dK1_Ir125,2)]);
hold on 
h = plot(x, a_Ir125*exp(-b_Ir125*x)+c_Ir125);  
h.Color = [0.9290 0.6940 0.1250]
xlim([8 800]); ylim([0.42 0.5])
% set(gca,'box','off')
set(gca,'LineWidth',1.7,'FontSize',12,'Box', 'on'); 
xlabel('Spin-Echo Time (ps)')
ylabel('Polarization')

subplot(2,3,4)
plot(dKIr125,alpha_Ir125,'bo','MarkerSize',6,'color',[.4 .4 .4]);
hold on
errorbar(dKIr125,alpha_Ir125,dalpha_Ir125,'LineStyle','None','color',[.4 .4 .4]);
ylim([-0.001 0.015]); xlim([0 3.2]);
hold on

%Fit alphas to the CE model 
[deltaks_Ir125, alphas_Ir125] = prepareCurveData(dKIr125,alpha_Ir125);
ft = fittype( '2/(t)*((p1*sin(x*2.72./2*cos(0)).^2)+(p1*sin(x*2.72./2*cos(60)).^2)+(p1*sin(x*2.72./2*cos(120)).^2)+(p1*sin(x*2.72./2*cos(180)).^2)+(p1*sin(x*2.72./2*cos(240)).^2)+(p1*sin(x*2.72./2*cos(300)).^2)+((1-p1)*sin(x*4.71./2*cos(30)).^2)+((1-p1)*sin(x*4.71./2*cos(90)).^2)+((1-p1)*sin(x*4.71./2*cos(150)).^2)+((1-p1)*sin(x*4.71./2*cos(210)).^2)+((1-p1)*sin(x*4.71./2*cos(270)).^2)+((1-p1)*sin(x*4.71./2*cos(330)).^2))', 'independent', 'x', 'dependent', 'y', 'coefficients',{'t','p1'} );                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 fittype( '2/(t)*((p1*sin(x*2.72./2*cos(30)).^2)+(p1*sin(x*2.72./2*cos(90)).^2)+(p1*sin(x*2.72./2*cos(150)).^2)+(p2*sin(x*4.71./2).^2)+(p2*sin(x*4.71./2*cos(60)).^2)+(p2*sin(x*4.71./2*cos(120)).^2))', 'independent', 'x', 'dependent', 'y', 'coefficients',{'t','p1','p2'} );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [65 0 ];
opts.StartPoint = [650 0.2 ];
opts.Upper = [10000 1];
opts.MaxFunEvals = 1000;
opts.MaxIter = 1000;
opts.TolFun = 1e-08;
opts.Weights = 1./dalpha_Ir125;
[fitresult_sin, gof_sin] = fit( deltaks_Ir125, alphas_Ir125, ft, opts )
fsin =(2/(fitresult_sin.t)*(fitresult_sin.p1*((sin(x_dKs*2.72./2*cos(0)).^2)+(sin(x_dKs*2.72./2*cos(60)).^2)+(sin(x_dKs*2.72./2*cos(120)).^2)+(sin(x_dKs*2.72./2*cos(180)).^2)+(sin(x_dKs*2.72./2*cos(240)).^2)+(sin(x_dKs*2.72./2*cos(300)).^2))+((1-fitresult_sin.p1)*((sin(x_dKs*4.71./2*cos(30)).^2)+(sin(x_dKs*4.71./2*cos(90)).^2)+(sin(x_dKs*4.71./2*cos(150)).^2)+(sin(x_dKs*4.71./2*cos(210)).^2)+(sin(x_dKs*4.71./2*cos(270)).^2)+(sin(x_dKs*4.71./2*cos(330)).^2)))));
plot(x_dKs,fsin,'color',[0.4660 0.6740 0.1880])
set(gca,'LineWidth',1.7,'FontSize',12,'Box', 'on');
xlabel('$\Delta K$(\AA$^{-1}$)');  ylabel('$\alpha$ (ps$^{-1}$)');

subplot(2,3,2)
GrIrCurve = semilogx(tseGrIr125_dK1,PmagGrIr125_dK1,'o','color',[0 0.4470 0.7410]);
hold on 
h = plot(x, a_GrIr125*exp(-b_GrIr125*x)+c_GrIr125)
h.Color = [0 0.4470 0.7410]
xlim([8 800]); ylim([0.1 0.18]);
%set(gca,'box','off')
set(gca,'LineWidth',1.7,'FontSize',12,'Box', 'on'); 
xlabel('Spin-Echo Time (ps)')



subplot(2,3,5)

plot(dKGrIr125,alpha_GrIr125,'bo','MarkerSize',6,'color',[.4 .4 .4]);
hold on
errorbar(dKGrIr125,alpha_GrIr125,dalpha_GrIr125,'LineStyle','None','color',[.4 .4 .4]);
xlim([0 3.2]);ylim([0 0.0005]);
set(gca,'LineWidth',1.7,'FontSize',12,'Box', 'on');
xlabel('$\Delta K$(\AA$^{-1}$)');

%figure(2)
subplot(2,3,3)
GrNiCurve = semilogx(tseGrNi125_dK1,PmagGrNi125_dK1,'diamond','color',[0.8500 0.3250 0.0980]);
hold on 
h = plot(x, Amp_GrNi125(index_GrNi125)*exp(-alpha_GrNi125(index_GrNi125)*x)+offset_GrNi125(index_GrNi125));
h.Color = [[0.8500 0.3250 0.0980]];
xlim([8 800])
ylim([0.2 0.28])
%set(gca,'box','off')
set(gca,'LineWidth',1.7,'FontSize',12,'Box', 'on'); 
%set(gcf,'Position',[0 0 1000 300])
xlabel('Spin-Echo Time (ps)')

subplot(2,3,6) %%Values copied from Tamtogl et al 
plot(dK_GrNi125,alpha_GrNi125,'bo','MarkerSize',6,'color',[.4 .4 .4]);
hold on
errorbar(dK_GrNi125,alpha_GrNi125,dalpha_GrNi125,'LineStyle','None','color',[.4 .4 .4]);
xlim([0 3.2]);
set(gca,'LineWidth',1.7,'FontSize',12,'Box', 'on');
xlabel('$\Delta K$(\AA$^{-1}$)');


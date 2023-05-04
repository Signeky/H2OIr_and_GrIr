

clear all; close all; 
addpath Dynamics\GrIr_and_Ir


%Load and fit data at 125 K 
fileListIr125 = {'dy019476.mat','dy019477.mat','dy019478.mat','dy019479.mat','dy019480.mat','dy019481.mat','dy019482.mat','dy019483.mat','dy019484.mat','dy019485.mat','dy019486.mat','dy019487.mat','dy019488.mat','dy019489'};  
for i=1:length(fileListIr125)
    
    % load in file from the list above
    load(fileListIr125{i})
    dKIr125(i) = meas.dK; 

    TIr125 = str2num(meas.endStatus.tSample);
    tseIr125 = meas.setime;
    PmagIr125 = meas.mean.Pmag;
    PimIr125 = meas.mean.Pimag;
    
    [xData, yData] = prepareCurveData(tseIr125(2:end),PmagIr125(2:end));
    ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.StartPoint = [0.4 0.01 min(PmagIr125)];
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
    
    %Store data and fit paramters at dK = 0.5 and dK = 0.7 
    if (dKIr125(i)<0.51) && (dKIr125(i)>0.49)
        alpha_Ir125_05 = alpha_Ir125(i); 
        dalpha_Ir125_05 = dalpha_Ir125(i);
        dK_Ir125_05 = dKIr125(i);
    elseif (dKIr125(i) < 0.71) && (dKIr125(i) > 0.69);    
        alpha_Ir125_07 = alpha_Ir125(i); 
        dalpha_Ir125_07 = dalpha_Ir125(i);
        dK_Ir125_07 = dKIr125(i);
    end 
  

   
end

%Load and fit data at 135 K 
fileListIr135 = {'dy019490.mat','dy019491.mat','dy019492.mat','dy019493.mat','dy019494.mat','dy019495.mat','dy019496.mat','dy019497.mat','dy019498.mat','dy019499.mat','dy019500.mat','dy019501.mat','dy019502.mat','dy019503'};  
for i=1:length(fileListIr135)
    
    % load in file from the list above
    load(fileListIr135{i})
    dKIr135(i) = meas.dK; 

    TIr135 = str2num(meas.endStatus.tSample);
    tseIr135 = meas.setime;
    PmagIr135 = meas.mean.Pmag;
    PimIr135 = meas.mean.Pimag;
    
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
    
    %Store data and fit paramters at dK = 0.5 and dK = 0.7
    if (dKIr135(i)<0.51) && (dKIr135(i)>0.49)
        alpha_Ir135_05 = alpha_Ir135(i); 
        dalpha_Ir135_05 = dalpha_Ir135(i);
        dK_Ir135_05 = dKIr135(i);
    elseif (dKIr135(i) < 0.71) && (dKIr135(i) > 0.69);    
        alpha_Ir135_07 = alpha_Ir135(i); 
        dalpha_Ir135_07 = dalpha_Ir135(i);
        dK_Ir135_07 = dKIr135(i);
    end 
  

   
end


%% Arrhenius plots and calculations of the activation energy 
%Data at 125 K 
dK_Ir125 = [dK_Ir125_05,dK_Ir125_07];
xc = 125; 
yc = [alpha_Ir125_05, alpha_Ir125_07].*1e3;
yc_error = [dalpha_Ir125_05, dalpha_Ir125_07].*1e3;

%Data at 135 K
dK_Ir135 = [dK_Ir135_05,dK_Ir135_07];
xh = 135; 
yh = [alpha_Ir135_05, alpha_Ir135_07].*1e3;
yh_error = [dalpha_Ir135_05, dalpha_Ir135_07].*1e3;

%Constants
BoltzmannConstant=1.3806504E-23;

%Prepare curve for dK = 0.5
x_05 = [1/xh, 1/xc];
y_05 = [log(yh(1)), log(yc(1))];
y_05_error = [yh_error(1)/yh(1), yc_error(1)/yc(1)]; %Error through error propagation 

slope_05 = (y_05(2)-y_05(1))/(x_05(2) - x_05(1));
slope_error_05 = sqrt((y_05_error(1)/(x_05(2)-x_05(1)))^2+(y_05_error(2)/(x_05(2)-x_05(1)))^2);
Ea_05 = abs(slope_05*BoltzmannConstant*6.24150974E21);
Ea_error_05 = abs(slope_error_05*BoltzmannConstant*6.24150974E21);
intercept_05 = y_05(1)-(slope_05*x_05(1));
alpha0_05 = exp(intercept_05)*1e-3; %in ps-1

%alpha_0  = alpha* e^(kx), find error through error proporgation. Choose alpha-error at 135 K?
alpha0_error_05_135 = sqrt((yh_error(1)*exp(Ea_05/(6.24150974E21*BoltzmannConstant*135)))^2+(Ea_error_05*yh(2)*exp(Ea_05/(6.24150974E21*BoltzmannConstant*135))/(BoltzmannConstant*6.24150974E21*135))^2)*1e-3;
alpha0_error_05_125 = sqrt((yc_error(1)*exp(Ea_05/(6.24150974E21*BoltzmannConstant*125)))^2+(Ea_error_05*yc(2)*exp(Ea_05/(6.24150974E21*BoltzmannConstant*125))/(BoltzmannConstant*6.24150974E21*125))^2)*1e-3;

fprintf(['At dK = 0.5, activation energy Ea = (%2.0f ' char(177) ' %2.0f) meV is found'], Ea_05, Ea_error_05)
fprintf([' with an exponential prefactor alpha_0 = (%2.0f ' char(177) ' %2.0f) ps^(-1)'], alpha0_05, alpha0_error_05_135)

%Prepare curve for dK = 0.7
x_07 = [1/xh, 1/xc];
y_07 = [log(yh(2)), log(yc(2))];
y_07_error = [yh_error(2)/yh(2), yc_error(2)/yc(2)]; %Error through error propagation 

slope_07 = (y_07(2)-y_07(1))/(x_07(2) - x_07(1));
slope_error_07 = sqrt((y_07_error(1)/(x_07(2)-x_07(1)))^2+(y_07_error(2)/(x_07(2)-x_07(1)))^2);
Ea_07 = abs(slope_07*BoltzmannConstant*6.24150974E21);
Ea_error_07 = abs(slope_error_07*BoltzmannConstant*6.24150974E21);

intercept_07 = y_07(1)-(slope_07*x_07(1));
alpha0_07 = exp(intercept_07)*1e-3; %in ps-1

%alpha_0  = alpha* e^(-kx), choose error at 135 K?
alpha0_error_07_135 = sqrt((yh_error(2)*exp(Ea_07/(6.24150974E21*BoltzmannConstant*135)))^2+(Ea_error_07*yh(2)*exp(Ea_07/(6.24150974E21*BoltzmannConstant*135))/(BoltzmannConstant*6.24150974E21*135))^2)*1e-3;
alpha0_error_07_125 = sqrt((yc_error(2)*exp(Ea_07/(6.24150974E21*BoltzmannConstant*125)))^2+(Ea_error_07*yc(2)*exp(Ea_07/(6.24150974E21*BoltzmannConstant*125))/(BoltzmannConstant*6.24150974E21*125))^2)*1e-3;

fprintf(['\nAt dK = 0.7, activation energy Ea = (%2.0f ' char(177) '%4.0f) meV is found'], Ea_07, Ea_error_07)
fprintf([' with an exponential prefactor alpha_0 = (%2.0f ' char(177) ' %2.0f) ps^(-1)'], alpha0_07, alpha0_error_07_135)

% Plot figure 
xslope=linspace(0.95*(x_07(1)),1.05*x_07(2),100);
figure=figure('units','centimeters','position',[5,3,14,10],'color','white','DefaultLineLineWidth',1.5);
set(figure,'DefaultLineLineWidth',2); set(gca,'LineWidth',1.5,'FontSize',14); box(gca,'on'); hold(gca,'all'); grid off;
hold on
plot(x_05,y_05, 'o','color',[0 0.4470 0.7410],'DisplayName', '\Delta K = 0.5');
errorbar(x_05,y_05,y_05_error,'color',[0 0.4470 0.7410],'LineStyle','None','HandleVisibility','off');
plot(xslope,slope_05*xslope+intercept_05, 'color',[0 0.4470 0.7410],'HandleVisibility','off');
plot(x_07,y_07, 'o','color',[0.4940 0.1840 0.5560],'DisplayName', '\Delta K = 0.7');
errorbar(x_07,y_07,y_07_error,'color',[0.4940 0.1840 0.5560],'LineStyle','None','HandleVisibility','off');
plot(xslope,slope_07*xslope+intercept_07, 'color',[0.4940 0.1840 0.5560],'HandleVisibility','off');
legend
xlim([0.0073 0.0081]); %ylim([-6.8 -4.8]);

xlabel('1/T (K^{-1})')
ylabel('ln(\alpha) (ns^{-1})')


%% Alternatively, the activation energy can be calculated using alpha_0 = (5 +- 1) ps^(-1) found for GrNi(111) 
%Error proporgation formula, Ea_error = sqrt(((sigmaalpha/alpha)6.24150974E21*BoltzmannConstant*125)^2+((sigmaalpha/alpha)*6.24150974E21*BoltzmannConstant*125)^2)

alpha0_GrNi = 5000; %in ns^(-1)
alpha0_GrNi_error = 1000; %in ns^(-1)

%At 0.5
slope_05_alternative = (y_05(2)-log(alpha0_GrNi))/x_05(2);
Ea_05_alternative = abs(slope_05_alternative*BoltzmannConstant*6.24150974E21);
Ea_05_alternative_error = sqrt(((yc_error(1)/yc(1))*6.24150974E21*BoltzmannConstant*125)^2+((alpha0_GrNi_error/alpha0_GrNi)*6.24150974E21*BoltzmannConstant*125)^2);

%At 0.7
slope_07_alternative = (y_07(2)-log(alpha0_GrNi))/x_07(2);
Ea_07_alternative = abs(slope_07_alternative*BoltzmannConstant*6.24150974E21);
Ea_07_alternative_error = sqrt(((yc_error(2)/yc(2))*6.24150974E21*BoltzmannConstant*125)^2+((alpha0_GrNi_error/alpha0_GrNi)*6.24150974E21*BoltzmannConstant*125)^2);

fprintf(['\nAltneratively, using the exponential prefactor of 5 ps(^-1) found for GrNi, we get'])
fprintf(['Ea = (%2.0f ' char(177) '%2.0f) at dK = 0.5 Å^(-1) and \n'],Ea_05_alternative,Ea_05_alternative_error )
fprintf(['Ea = (%2.0f ' char(177) '%2.0f) at dK = 0.7 Å^(-1)\n'],Ea_07_alternative,Ea_07_alternative_error )


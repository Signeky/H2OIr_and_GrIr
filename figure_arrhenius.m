

clear all; close all; 
addpath Dynamics\Ir_and_GrIr

%% First load and fit data 
%At 125 K 
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
    
    %Store data and fit paramters at dK = 0.7 

    if (dKIr125(i) < 0.71) && (dKIr125(i) > 0.69);    
        alpha_Ir125_07 = alpha_Ir125(i); 
        dalpha_Ir125_07 = dalpha_Ir125(i);
        dK_Ir125 = dKIr125(i);

    end 
  

   
end

%At 135 K 
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
    
    %Store data and fit paramters at dK = 0.7
    if (dKIr135(i) < 0.71) && (dKIr135(i) > 0.69);    
        alpha_Ir135_07 = alpha_Ir135(i); 
        dalpha_Ir135_07 = dalpha_Ir135(i);
        dK_Ir135 = dKIr135(i);

    end 
  

   
end


%% Arrhenius plots and calculations of the activation energy 
%Data at 125 K 
xc = 125; 
yc = alpha_Ir125_07*1e3;
yc_error = dalpha_Ir125_07*1e3;

%Data at 135 K
xh = 135; 
yh = alpha_Ir135_07*1e3;
yh_error = dalpha_Ir135_07*1e3;

%Constants
BoltzmannConstant=1.3806504E-23;


%Prepare curve for dK = 0.7
x = [1/xh, 1/xc];
y = [log(yh), log(yc)];
y_error = [yh_error/yh, yc_error/yc]; %Error through error propagation 

%Calculate slope
slope = (y(2)-y(1))/(x(2) - x(1));
slope_error = sqrt((y_error(1)/(x(2)-x(1)))^2+(y_error(2)/(x(2)-x(1)))^2);
%Calculate activation energy 
Ea = abs(slope*BoltzmannConstant*6.24150974E21);
Ea_error = abs(slope_error*BoltzmannConstant*6.24150974E21);

intercept = y(1)-(slope*x(1));
alpha0 = exp(intercept)*1e-3; %in ps-1

%alpha_0  = alpha* e^(-kx)
alpha0_error_135 = sqrt((yh_error*exp(Ea/(6.24150974E21*BoltzmannConstant*135)))^2+(Ea_error*yh*exp(Ea/(6.24150974E21*BoltzmannConstant*135))/(BoltzmannConstant*6.24150974E21*135))^2)*1e-3;

%Print output
fprintf(['\nAt dK = 0.7, activation energy Ea = (%2.0f ' char(177) '%4.0f) meV is found'], Ea, Ea_error)
fprintf([' with an exponential prefactor alpha_0 = (%2.0f ' char(177) ' %2.0f) ps^(-1)'], alpha0, alpha0_error_135)

% Plot figure 
xslope=linspace(0.95*(x(1)),1.05*x(2),100);
figure=figure('units','centimeters','position',[5,3,14,10],'color','white','DefaultLineLineWidth',1.5);
set(figure,'DefaultLineLineWidth',2); set(gca,'LineWidth',1.5,'FontSize',14); box(gca,'on'); hold(gca,'all'); grid off;
hold on

plot(x,y, 'o','color',[0.4940 0.1840 0.5560],'DisplayName', '\Delta K = 0.7');
errorbar(x,y,y_error,'color',[0.4940 0.1840 0.5560],'LineStyle','None','HandleVisibility','off');
plot(xslope,slope*xslope+intercept, 'color',[0.4940 0.1840 0.5560],'HandleVisibility','off');
%legend
xlim([0.0073 0.0081]); ylim([1.1 2.25]);

xlabel('1/T (K^{-1})')
ylabel('ln(\alpha) (ns^{-1})')


%% Alternatively, the activation energy can be calculated using alpha_0 = (5 +- 1) ps^(-1) found for GrNi(111) 

alpha0_GrNi = 5000; %in ns^(-1)
alpha0_GrNi_error = 1000; %in ns^(-1)


slope_alternative = (y(2)-log(alpha0_GrNi))/x(2);
Ea_alternative = abs(slope_alternative*BoltzmannConstant*6.24150974E21);
Ea_alternative_error = sqrt(((yc_error/yc)*6.24150974E21*BoltzmannConstant*125)^2+((alpha0_GrNi_error/alpha0_GrNi)*6.24150974E21*BoltzmannConstant*125)^2);

fprintf(['\nAltneratively, using the exponential prefactor of 5 ps(^-1) found for GrNi, we get'])
fprintf([' Ea = (%2.0f ' char(177) '%2.0f) at dK = 0.7 Å^(-1)\n'],Ea_alternative,Ea_alternative_error )


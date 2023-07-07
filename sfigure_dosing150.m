% use the following to determine the pre-time
%
% clf; load('sm006786.mat'); semilogy(meas.t,meas.samp_pressure_e,meas.t,meas.counts) 
% 
clear all; clf


%File info
dataFile = 'Sm006788.mat';
temperature = '150 K';
pre_t = 45;


%Define figure
figure1=figure('units','centimeters','position',[5,3,18,18],'color','white','DefaultTextInterpreter','LaTex');
subplot(2,1,1);
set(gca,'LineWidth',1.5,'FontSize',14,'Layer','top'); box(gca,'on'); hold(gca,'all'); ax1=gca;
subplot(2,1,2);
set(gca,'LineWidth',1.5,'FontSize',14,'Layer','top','YScale','Log'); box(gca,'on'); hold(gca,'all'); ax2=gca;



%Get data
load(dataFile)
pressure = meas.samp_pressure_e;
pressure = (pressure(1:end-1)+pressure(2:end))/2;

time = (meas.t(1:end-1)+meas.t(2:end))/2;

signal = (meas.counts(1:end-1)+meas.counts(2:end))/2;
signal = signal./max(signal);

axes(ax1); 
plot(time,signal,'LineWidth',2.0)
xlabel('Time (s)'); ylabel('Normalised counts (arb units)');

axes(ax2);
semilogy(time,pressure,'LineWidth',2.0)
xlabel('Time (s)'); ylabel('Extractor pressure (mbar)');
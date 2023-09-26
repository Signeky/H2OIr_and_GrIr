
clear all; clf

fileList = {'Diffraction/As010351.mat','Diffraction/As010354'};
names = {'GrIr(111)','H_2O GrIr(111)'};

colourList = 'rbgcmk';

figure=figure('units','centimeters','position',[4,2,16,14],'color','white','DefaultLineLineWidth',2,'DefaultTextFontSize',18);

addpath('Diffraction', 'Calibration')


for i=1:length(fileList)
    load(fileList{i})
    signal = (meas.scan.counts(1:end-1)+meas.scan.counts(2:end))/2;
    signal = signal./max(signal);
    name = names{i};
    semilogy(meas.scan.angles,meas.scan.counts,colourList(i),'displayname', [num2str(name)]) %Plot angles as relative to the specular angle
    
    hold on 
end

xlabel('Angle of incidence, \theta_i(\circ)')
ylabel('Scattered signal (arb units)')
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
xlim([21 53])
ylim([4*10^(-11) 1.5*10^(-8)])

legend
grid
hold off



for i=1:length(fileList)
    load(fileList{i})
    signal = (meas.scan.counts(1:end-1)+meas.scan.counts(2:end))/2;
    signal = signal./max(signal);
    name = names{i};
    semilogy(meas.scan.angles-meas.ga,meas.scan.counts,colourList(i),'displayname', [num2str(name)]) %Plot angles as relative to the specular angle
    
    hold on 
end

xlabel('\Delta\theta_i relative to specular (\circ)')
ylabel('Scattered signal (arb units)')
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
xlim([21-meas.ga 53-meas.ga])
ylim([4*10^(-11) 1.5*10^(-8)])

legend
grid
hold off

%% Plot with Delta K rather than angle

% Use the calibration from 08/12/2021
load('beam_20211208_1219.mat')
E0 = beam.E0; % meV

E0 = E0*1.603e-19/1e3; % J
h = 6.626e-34; % Planck constant, Js
m = 3.0160293*1.66054e-27; % Helium 3 rest mass, kg
k = (2*pi*sqrt(2*E0*m)/h); % m^-1, incident wavector
k = k/1e10; % A^-1
disp(['k, incident wavevector, = ' num2str(k) ' A^-1'])

load(fileList{1})
IrGr_gamma = meas.scan.angles-meas.ga;
IrGr_counts = meas.scan.counts;
IrGr_DK = -2*k*cosd(22.2)*sind(IrGr_gamma);
load(fileList{2})
IrGrH20_gamma = meas.scan.angles-meas.ga;
IrGrH20_counts = meas.scan.counts;
IrGrH20_DK = -2*k*cosd(22.2)*sind(IrGrH20_gamma);

semilogy(-IrGr_DK, IrGr_counts, colourList(1), 'displayname', names{1})
hold on
semilogy(-IrGrH20_DK, IrGrH20_counts, colourList(2), 'displayname', names{2})
hold off
xlabel('\Delta K (Å^{-1})')
ylabel('Scattered signal (arb units)')
legend
set(gca,'YScale','log','LineWidth',1.5,'FontSize',18); set(gca,'ticklength',1.5*get(gca,'ticklength'))
grid on
xlim([0.11 3.4])
ylim([4*10^(-11) 1.5*10^(-8)])

% Find the diffraction peak to measure the lattice constant
ind = IrGr_DK < -2.8 & IrGr_DK > -3.16;
y = flip(IrGr_counts(ind));
x = flip(IrGr_DK(ind));
pks = findpeaks(y, x);
G = abs(x(y == pks));
disp(['Diffraction peak for graphene lattice spacing occurs at ' num2str(G) ' Å^-1'])
a = (2*pi/G)/cosd(30); % Account for measuring along the Gamma-M direction
disp(['Measured lattice constant of GrIr(111) = ' num2str(a) ' A'])





%please read the readme file
% Code written by Kamran Keykhosravi oct 2020

close all;clear all;clc;
rand_seed =1;
%% generating ris profile
config = GetConfig();
rng(100)
ris_phases = exp(1j*2*pi*rand(config.Mc^2,config.T));
%% EStimation and CRB calculation
estimator = Estimator(config,ris_phases,rand_seed);
estimator = estimator.estimate; % calculating peb
estimator = estimator.PEBcalc; % calculating crb bounds
% calculating scaling of peb with number of ris elements
config2=config; 
config2.PointNum=4;
config2.r_vec = [20,30,40,50];
config2.xyz = [-config2.r_vec/sqrt(2);config2.r_vec/sqrt(2);-10*ones(1,length(config2.r_vec))];
config2.McVec = floor(10.^[1:.2:3]);
PEB = zeros(config2.PointNum,length(config2.McVec));
for imc = 1:length(config2.McVec)
    fprintf('Total %.2f %% done! \n', (imc-1)/length(config2.McVec)*100);
    config2.Mc = config2.McVec(imc);
    estimator2 = Estimator(config2,exp(1j*2*pi*rand(config2.Mc^2,config2.T)),rand_seed);
    estimator2 = estimator2.PEBcalc; % calculating crb bounds
    PEB(:,imc)= estimator2.PEB.';
end

%% plotting
figure %position error bound and \Delta_t
plot(estimator.config.r_vec,estimator.PEB)
hold on
plot(estimator.config.r_vec,sqrt(mean(estimator.Error_Squared.PEB,1)),'x')
plot(estimator.config.r_vec,config.c*estimator.getCrb(estimator.FIMpo,4),'--')
plot(estimator.config.r_vec,config.c*sqrt(mean(estimator.Error_Squared.DT,1)),'+')
set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');
xlabel('r (m)');ylabel('Error(m)');legend('CRB Pos','Est Pos','CRB \Delta_t','Est. \Delta_t');

%

figure %\tau_b and \tau_r
plot(estimator.config.r_vec,config.c*estimator.getCrb(estimator.FIMch,1))
hold on
plot(estimator.config.r_vec,config.c*sqrt(mean(estimator.Error_Squared.tau_b,1)),'x')
plot(estimator.config.r_vec,config.c*estimator.getCrb(estimator.FIMch,2),'--')
plot(estimator.config.r_vec,config.c*sqrt(mean(estimator.Error_Squared.tau_r,1)),'+')
set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');
xlabel('r (m)');ylabel('Error(m)');legend('CRB \tau_b','Est \tau_b','CRB \tau_r','Est. \tau_r');

%

figure %\phi_az and \phi_el
plot(estimator.config.r_vec,config.c*estimator.getCrb(estimator.FIMch,3))
hold on
plot(estimator.config.r_vec,config.c*sqrt(mean(estimator.Error_Squared.phi_az,1)),'x')
plot(estimator.config.r_vec,config.c*estimator.getCrb(estimator.FIMch,4),'--')
plot(estimator.config.r_vec,config.c*sqrt(mean(estimator.Error_Squared.phi_el,1)),'+')
set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');
xlabel('r (m)');ylabel('Error(m)');legend('CRB \phi_{az}','Est \phi_{az}','CRB \phi_{el}','Est. \phi_{el}');

%

figure
plot(config2.McVec.^2,PEB)
set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');
xlabel('Number of RIS elements');ylabel('PEB(m)');legend('r=20','r=30','r=40','r=50');


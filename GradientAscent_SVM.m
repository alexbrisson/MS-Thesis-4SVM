%% Gradient Ascent optimization - Main
% This script connects to OMA software and runs the optimization of 4-SVM
% locally. 

% Alex Brisson
% 3/24
%% Transmitter Parameters
if ~exist('txparams')
txparams = 'are loaded';
tx.dacClk = 92e9; % hw specific, some flexbility if internally generated
tx.desBaud = 4.6e9; % desired Baud Rate
tx.desSapSym = tx.dacClk/tx.desBaud; % desired samples per symbol
tx.M = 4; % SVM Constellation Cardinality
tx.NominalJonesVectors = txJonesVectors(1);
tx.NominalStokesVectors = jones2stokes_ab(tx.NominalJonesVectors)';
tx.deB = csvread('prss_4_5.csv'); % load in de Bruijn Sequence

% RRC Filter parameters if RRC is used
tx.rrc_rolloff = 1;
tx.rrc_timespan = 10; 
% Generate AWG Waveforms
[tx.sqWaveforms,tx.vdack,tx.rrcWaveforms,tx.voltSeq] = svm_AWGWaveforms(tx.NominalJonesVectors, ...
    tx.desSapSym,tx.deB,tx.rrc_rolloff,tx.rrc_timespan);

% Plot Waveforms
if 0 
    figure()
    plot(tx.sqWaveforms(:,1),'color','k','LineWidth',4);hold on
    plot(tx.sqWaveforms(:,2),'color','m','LineWidth',4)
    plot(tx.sqWaveforms(:,3),'color','g','LineWidth',4)
    plot(tx.sqWaveforms(:,4),'color','b','LineWidth',4)
    ylabel('V\pi');xlabel('sample');legend('xi','xq','yi','yq')
    xlim([0,700]);ylim([-1.1,1.1]);
    grid on;hold off
end
end

% Optimization plot, real-time updates
figure(2)
plt1 = plot(0,'o','color','k'); hold on
plt2 = plot(0,'o','color','r');
title('results');xlabel('Iteration Index');grid on
legend('\xi/\xi_{max}','EVM')

%% Establish Connection to OMA and AWG

sysparams.objscurr = instrfind;
for i=1:length(sysparams.objscurr)
    if strcmp('open',sysparams.objscurr.status(i))
        fclose(sysparams.objscurr(i));
    end
end
oma = WCFServiceOM4006Basic; % obj for oma control functions
awg = visa('AGILENT', 'TCPIP0::localhost::inst0::INSTR');
fopen(awg)

% set initial run iteration
runiter = 1;

%% Main Optimization Loop

while true 

% -- check if oui is running, if so, wait. 
sysparams.omaStateProcessing =1; % set state to running true
while sysparams.omaStateProcessing ==1
pause(0.1)
sysparams.omaIsProcessingbool = IsProcessing(oma);
sysparams.omaStateProcessing = strcmp('true',sysparams.omaIsProcessingbool);
end

Single(oma) % take single aquisition from OUI (aquisition length is based on record length set in OUI)

% -- Wait for single aquisition to finish processing
sysparams.omaStateProcessing =1; % set state to running true
while sysparams.omaStateProcessing ==1
pause(0.1) % reduce frequency of checks, don't clog up oma
sysparams.omaIsProcessingbool = IsProcessing(oma);
sysparams.omaStateProcessing = strcmp('true',sysparams.omaIsProcessingbool);
end

while ~exist('evm.mat')
    pause(0.1)
end

load('stokesn.mat'); % load in normalized cluster averaged Stokes vectors
load('evm.mat') % load in EVM measured

% delete .mat files
delete('stokesn.mat'); 
delete('evm.mat')

% turn AWG output off
fprintf(awg,':abor');

% set k-1-th iteration to k-th current
if runiter == 1
    jonesnk_1 = tx.NominalJonesVectors; % Jones Vectors
    vdack_1 = tx.vdack; % Gamma Values (normalized voltages)
else
    jonesnk_1 = jonesnk;
    vdack_1 = vdack;
end

% Gradient Ascent Function to Optimize Voltages--- Outputs next iteration
% of gamma values with cost function value xi.
step = 1e-3; % step size for optimization
[xi,jonesnk,stokesnk,vdack] = voltOpt(mlc2,jonesnk_1,vdack_1,step);
disp(num2str(vdack));

%% Generate new set of AWG waveforms using k-th gamma values

szprss = size(tx.deB);
sapsym = 20;
vseq = zeros([max(szprss),4]);
for i=1:max(szprss)
    vseq(i,:) = vdack(tx.deB(i),:);
end

% Interpolate voltages to samples/symbol
for i = 1:length(vseq)
vnmseq_i=repmat(vseq(i,:),sapsym,1);
if i ==1
    vnmseq = vnmseq_i;
else
    vnmseq = [vnmseq;vnmseq_i];
end
end

% write waveforms to current working folder
csvwrite('v1.csv',vnmseq(:,1));
csvwrite('v2.csv',vnmseq(:,2));
csvwrite('v3.csv',vnmseq(:,3));
csvwrite('v4.csv',vnmseq(:,4));

%% Update optimization plot // write new waveforms to AWG

if runiter == 1
    pltxi = xi;
    pltevm = evm;
else
    pltxi = [pltxi,xi];
    pltevm = [pltevm,evm];
end

set(plt1,'YData',pltxi./9.48) % xi normalized to maximum
set(plt2,'YData',pltevm)

% ensure new waveforms are saved to workspace
while ~exist('v4.csv')
    pause(0.1)
end

% Upload new waveforms to AWG channels
fprintf(awg,strcat([':trac1:imp 1, "',pwd,'\v1.csv",CSV,IONL,OFF,ALEN']));
fprintf(awg,strcat([':trac2:imp 1, "',pwd,'\v2.csv",CSV,IONL,OFF,ALEN']));
fprintf(awg,strcat([':trac3:imp 1, "',pwd,'\v3.csv",CSV,IONL,OFF,ALEN']));
fprintf(awg,strcat([':trac4:imp 1, "',pwd,'\v4.csv",CSV,IONL,OFF,ALEN']));

% Delete the waveforms from working folder
pause(0.1)
fclose('all')
delete('v1.csv');
delete('v2.csv');
delete('v3.csv');
delete('v4.csv');

% ensure new waveforms are deleted from workspace
while exist('v4.csv')/2
    pause(0.1)
end

% turn on AWG increment iteration
fprintf(awg,':init:imm');
pause(0.2)
runiter = runiter+1;
end
% Core processing code for SVM
% This script runs DSP algorithms on the OMA software

% Alex B.
% 03/24
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

%% Receiver Parameters

rx.recordLength = BlockSize; % should be long compared to ADC rate*PRSS length to allow for sufficient delay resolution
rx.adcClock = 100e9;
rx.timeSpan = (1/(rx.adcClock))*rx.recordLength;
rx.Nsymbols = (rx.timeSpan)*tx.desBaud;
rx.Sps = ceil(rx.adcClock/tx.desBaud); % integer number of samples per symbol
rx.interpolatedRecordLength = rx.Sps*rx.Nsymbols;

%% Set up alerts
if ~exist('ZonesAlreadyCalledFrom','var');ZonesAlreadyCalledFrom = [];end;
AllAlerts = [];
AllAlerts = [AllAlerts,struct('Code',021,'Msg','Cannot execute 2nd phase estimate because one or more tribs is not synchronised.')];
FunctionName = 'CoreProcessing';
[AllAlerts.FunctionName] = deal(FunctionName);
CurrZone = '';
if isempty(ZonesAlreadyCalledFrom)
   Alerts = RegisterAlerts(AllAlerts,Alerts);
   ZonesAlreadyCalledFrom = {[CurrZone,'x']};
elseif ~any(strcmp(ZonesAlreadyCalledFrom,[CurrZone,'x']))
   Alerts = RegisterAlerts(AllAlerts,Alerts);
   ZonesAlreadyCalledFrom = [ZonesAlreadyCalledFrom;[CurrZone,'x']];
end

% Initialise block processing boundary values
if FirstBlock
   EstimateClockBoundVals = [];
   ClockRetimeBoundVals1 = [];
   EstimateSOPBoundVals1 = [];
   EstimatePhaseBoundVals1 = [];
   ApplyPhaseBoundVals1 = [];
   QDecThBoundVals1 = [];
   QDecThBoundVals2 = [];
   ClockRetimeBoundVals2 = [];
   ApplyPhaseBoundVals2 = [];
   TrimTimeGridBoundVals1 = [];
   AlignTribsBoundVals1 = [];
   GenPatternBoundVals9 = [];
   GenPatternBoundVals10 = [];
   EstimatePhaseBoundVals3 = [];
   ApplyPhaseBoundVals7 = [];
   ApplyPhaseBoundVals8 = [];
   TrimTimeGridBoundVals2 = [];
   GenPatternBoundVals13 = [];
   GenPatternBoundVals14 = [];
   DiffDetectionBoundVals6 = [];
   ApplyPhaseBoundVals3 = [];
   ApplyPhaseBoundVals4 = [];
   TrimTimeGridBoundVals3 = [];
end

%% Matched Filter
Alerts.CurrZone = 'Impairment compensation';
% Correct for ch dealy and hybrid imbalance before filtering. 
[pVirt,ClockRetimeBoundVals1,Alerts] = ClockRetime(Vblock,ChDelay,pHyb,Vblock(1),ClockRetimeBoundVals1,Alerts);
pVirt_f = CopyTime(pVirt);
% rx matched filter impulse response
rx.h=(1/round(tx.desSapSym))*ones(1,round(tx.desSapSym));
pVirt_f_x = circshift(upfirdn(pVirt.Values(1,:),rx.h),-10,2);
pVirt_f_y = circshift(upfirdn(pVirt.Values(2,:),rx.h),-10,2);
% filtered rx waveforms
pVirt_f.Values = [pVirt_f_x(1:length(pVirt.Values));pVirt_f_y(1:length(pVirt.Values))];
% new block structure
V = [struct('t0',pVirt_f.t0,'dt',pVirt_f.dt,'Values',real(pVirt_f.Values(1,:))), ...
         struct('t0',pVirt_f.t0,'dt',pVirt_f.dt,'Values',imag(pVirt_f.Values(1,:))), ...
         struct('t0',pVirt_f.t0,'dt',pVirt_f.dt,'Values',real(pVirt_f.Values(2,:))), ...
         struct('t0',pVirt_f.t0,'dt',pVirt_f.dt,'Values',imag(pVirt_f.Values(2,:)))];
% already compensated- use in clk recovery
ChDelay1 = [0,0,0,0];
pHyb1 = [1,1i,0,0;0,0,1,1i];

%% Clock recovery
Alerts.CurrZone = 'Clock recovery';
NonlinFunc = '';
[SymClock,EstimateClockBoundVals,Alerts] = ...
   EstimateClock(V,ChDelay1,pHyb1,FreqWindow,NonlinFunc,EstimateClockBoundVals,Alerts);
TraceClock.t0 = SymClock.t0;
TraceClock.dt = SymClock.dt/TracePtsPerSym;
[pSym,ClockRetimeBoundVals1,Alerts] = ...
   ClockRetime(V,ChDelay1,pHyb1,SymClock,ClockRetimeBoundVals1,Alerts);
pSym.t = pSym.t0:pSym.dt:(pSym.dt*(length(pSym.Values)-length(pSym.t0:pSym.dt:0)));

% In-phase-Quadrature scatter plot
if 0
    figure();
    scatter(real(pSym.Values(1,:)),imag(pSym.Values(1,:)));hold on
    scatter(real(pSym.Values(2,:)),imag(pSym.Values(2,:)),'.');
    title('xpol/ypol clk retime');legend('xpol','ypol');axis equal;hold off
end

%% Jones 2 Stokes
% Convert Jones to Stokes and Normalize 
[pSymSt,pSymNmSt] = CopyTime(pSym);
pSymSt.Values = Jones2Stokes(pSym.Values); % non-normalized stokes vectors
pSymNmSt.Values = StokesNorm(pSymSt.Values); % normalized stokes vectors
rx.NsymbolsActual = length(pSymNmSt.Values); % actual Num Symbols Received in record

% Plot Nm Stokes Vectors on the Poincare Sphere
if 0
    figure()
    [x,y,z]=sphere;
    surf(x,y,z,'FaceAlpha',.1,'EdgeColor','none');hold on
    plot3(pSymNmSt.Values(1,:),pSymNmSt.Values(2,:),pSymNmSt.Values(3,:),'.')
    axis equal
end

%% Symbol Sequence Resolution

% Resolve Tx/Rx Symbol Delay
% Transmitted voltage sequence, jones sequence, and voltage indices for symbol delay res. 
tx.repindex = ceil(rx.NsymbolsActual/length(tx.deB)); % repeat symbol sequence
tx.vindseq = repmat(tx.deB(:,1),[tx.repindex,1]);
[tx.vindseq_crot,rx.symDelay,rx.xcorrs,mlc] = ...
    svm_SymDelayRes_ab(pSymNmSt.Values,tx.vindseq,tx.M);
bitPRSS = tx.vindseq_crot(1:rx.NsymbolsActual);
txbits = transbitsGen_ab(bitPRSS,tx.M);
tx.stokesSeq = tx.NominalStokesVectors(bitPRSS,:);
tx.jonesSeq = tx.NominalJonesVectors(bitPRSS,:);

%% maximum likelihood used to align symbols
% Cluster Averages: mlc2

mlc2 = zeros(size(mlc));
for i1 = 1:tx.M
    tempInds = find(bitPRSS == i1);
    decnormInd=zeros([1,length(tempInds)]);
    for i2 = 1:length(tempInds)
    [~,decnormInd(i2)] = max(pSymNmSt.Values(:,tempInds(i2))'*mlc');
    end
    decnormInd_likely = mode(decnormInd);
    mlc2(i1,:) = mlc(decnormInd_likely,:);
end

%% Pseudo Inverse

RotM = (pSymNmSt.Values*tx.stokesSeq)/(tx.stokesSeq'*tx.stokesSeq);
mlc2_crot = (inv(RotM)*mlc2')';
pSymNmSt_crot.Values = inv(RotM)*pSymNmSt.Values;
for i = 1:tx.M
    mlc2_crot(i,:) = mlc2_crot(i,:)./(norm(mlc2_crot(i,:)));
end
for i = 1:length(pSymNmSt.Values)
    pSymNmSt_crot.Values(:,i) = pSymNmSt_crot.Values(:,i)./norm(pSymNmSt_crot.Values(:,i));
end

%% Compute BER
[rxBER.BER,rxBER.TotalBits,rxBER.NumErrs] = BEReval_ab(tx.NominalStokesVectors,pSymNmSt_crot.Values,tx.M,txbits);
rxBER
    
%% Plotting Visualizations//Statistics

% compute EVM
[evm,mlcf] = SVM_evm(pSymNmSt_crot.Values,tx.M);

if strcmp(poincarefig,'true')
figure(1)
[x,y,z]=sphere;

rsphere_norm = zeros([1,length(pSymSt.Values)]);
for i = 1:length(pSymSt.Values)
rsphere_norm(i) = norm(pSymSt.Values(:,i));
end

rsphere_norm = mean(rsphere_norm);

r = rsphere_norm;
xn = x*r;
yn = y*r;
zn = z*r;
if strcmp(poincareenable_stnorm,'true')
    surf(x,y,z,'FaceAlpha',.1,'EdgeColor','none');
else
    surf(xn,yn,zn,'FaceAlpha',.1,'EdgeColor','none');
end
hold on
% plots for the Stokes points
plt1=plot3(0,0,0,'.','MarkerSize',12,'color','b');
plt2=plot3(0,0,0,'.','MarkerSize',12,'color','b');

% plots for the mean stokes vectors
plt3=plot3(0,0,0,'o','color','k','LineWidth',10);
plt4=plot3(0,0,0,'o','color','r','LineWidth',10);
plt5=plot3(0,0,0,'o','color','m','LineWidth',10);
plt6=plot3(0,0,0,'o','color','g','LineWidth',10);

plt7=plot3(0,0,0,'x','color','k','LineWidth',10);
plt8=plot3(0,0,0,'x','color','r','LineWidth',10);
plt9=plot3(0,0,0,'x','color','m','LineWidth',10);
plt10=plot3(0,0,0,'x','color','g','LineWidth',10);
axis equal
hold off
end

mnplt = mlc2_crot;

% plot values
if strcmp(poincareenable_stnorm,'true')
set(plt1,'Xdata',pSymNmSt_crot.Values(1,:));
set(plt1,'Ydata',pSymNmSt_crot.Values(2,:));
set(plt1,'Zdata',pSymNmSt_crot.Values(3,:));

set(plt3,'Xdata',tx.NominalStokesVectors(1,1));
set(plt3,'Ydata',tx.NominalStokesVectors(1,2));
set(plt3,'Zdata',tx.NominalStokesVectors(1,3));

set(plt4,'Xdata',tx.NominalStokesVectors(2,1));
set(plt4,'Ydata',tx.NominalStokesVectors(2,2));
set(plt4,'Zdata',tx.NominalStokesVectors(2,3));

set(plt5,'Xdata',tx.NominalStokesVectors(3,1));
set(plt5,'Ydata',tx.NominalStokesVectors(3,2));
set(plt5,'Zdata',tx.NominalStokesVectors(3,3));

set(plt6,'Xdata',tx.NominalStokesVectors(4,1));
set(plt6,'Ydata',tx.NominalStokesVectors(4,2));
set(plt6,'Zdata',tx.NominalStokesVectors(4,3));

set(plt7,'Xdata',mnplt(1,1));
set(plt7,'Ydata',mnplt(1,2));
set(plt7,'Zdata',mnplt(1,3));

set(plt8,'Xdata',mnplt(2,1));
set(plt8,'Ydata',mnplt(2,2));
set(plt8,'Zdata',mnplt(2,3));

set(plt9,'Xdata',mnplt(3,1));
set(plt9,'Ydata',mnplt(3,2));
set(plt9,'Zdata',mnplt(3,3));

set(plt10,'Xdata',mnplt(4,1));
set(plt10,'Ydata',mnplt(4,2));
set(plt10,'Zdata',mnplt(4,3));
end

if strcmp(poincareenable_st,'true')
set(plt2,'Xdata',pSymSt.Values(1,:));
set(plt2,'Ydata',pSymSt.Values(2,:));
set(plt2,'Zdata',pSymSt.Values(3,:));
end





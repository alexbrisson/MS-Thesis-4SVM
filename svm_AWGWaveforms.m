function [vnmseq,P,vnmseqRRC,vseq] = svm_AWGWaveforms(jones,sapsym,prss,rolloff,span)
%% This Function Generates the Normalized to Vpi Voltage Sequence for M-ary SVM
% By Alex Brisson, Optical Communications Group
% 01/24

%% Inputs: 
% jones: Symbol Jones Vectors size Mx2 
% sapsym: samples per symbol 
% prss: desired pseudo random symbol sequence
% rolloff: rolloff factor for rrc pulses
% span: span of rrc pulses

% symbol jones vectors should be size Mx2
szjones = size(jones);
if szjones(2) ~= 2
    jones = jones.';
    szjones = size(jones);
end
szprss = size(prss);
beta = pi/2; % FM Index
dac_sc = 1; % maximum desired DAC scaling ie, [0,1] (use if a smaller dynamic range is needed)
jv_sc = sin((beta)*dac_sc);

%Max scale the Jones vector set; sign the symbol sequence for AMI coding.
for i = 1:szjones(1)
    jones_sc(i,:)=jones(i,:).*(jv_sc/max([abs(real(jones(i,1))),abs(imag(jones(i,1))), ...
        abs(real(jones(i,2))),abs(imag(jones(i,2)))]));
    symb_ind=find(i==prss);
end

% Solve for gamma values
for i = 1:szjones(1)
    Pxi(i,1) = asin(real(jones_sc(i,1)))/beta;
    Pxq(i,1) = asin(imag(jones_sc(i,1)))/beta;
    Pyi(i,1) = asin(real(jones_sc(i,2)))/beta;
    Pyq(i,1) = asin(imag(jones_sc(i,2)))/beta;
end
P = [Pxi,Pxq,Pyi,Pyq];

for i = 1:szprss(1)
    vseq(i,:) = P(prss(i,1),:);
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

% generate the RRC pulse and remove the time delay
rrc_taps = rcosdesign(rolloff,span,sapsym);
rrc_taps = rrc_taps./(max(rrc_taps)); % scale to unit amplitude
vnmseqRRC = upfirdn(vseq,rrc_taps,sapsym);
half_IR_delay = round((length(vnmseqRRC) - length(vnmseq))/2);
vnmseqRRC(1:half_IR_delay,:)=[];
vnmseqRRC(length(vnmseq)+1:end,:)=[];
end
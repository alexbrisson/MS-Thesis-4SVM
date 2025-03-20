function [voltGroups_crot,symDelay,dval,mlc] = svm_SymDelayRes_ab(normSyms,vindseq,M)

%% Symbol Delay Resolution Function for any M-ary SVM Signal
% Alex Brisson, Optical Communications Group
% 11/2023

% This function resolves symbol delay between transmitter and reciever
% assuming transmitted voltage sequence is known. 

%% Inputs
% normSyms  = Normalized Symbol Sequence (normalized stokes vectors) 3xN
% vindseq   = voltage sequence from PRSS
% M         = Symbol Cardinality of SVM Signal (2-pol) 

%% Function Execution

voltGroups = vindseq;
% create a matrix list of all permuations of symbol integers
% i.e., there will be M! permutations.
tx.symbPerms = perms(1:M);

% Separate out each normalized stoke vector island and set each symbol
% vector to the corresponding integer in the nth permutation row.
% crosscorrelate the nth permutated symbol integer sequence with the voltage
% sequence until a correlation is resolved. 

tx.stokesNm = normSyms;
tx.szst = size(tx.stokesNm);
if tx.szst(1) == 3
tx.stokesNm = tx.stokesNm';
end

%% Seperate out groups of points and find the averages, normalize

Z = linkage(normSyms','ward'); % linkage: Agglomerative hierarchical cluster tree for clustering
c = cluster(Z,'Maxclust',M); % clusters groups into M groups.

for i1 = 1:M
    tx.stokesStr = strcat('tx.stokes_g',num2str(i1));
    gInds = find(c==i1);
    eval([tx.stokesStr,'=transpose(normSyms(:,gInds))'])
end

% Compute average Stokes vector for each cluster
for i1 = 1:M
    tx.stokesStr = strcat('tx.stokes_g',num2str(i1));
    mlc(i1,1:3) = [mean(eval(strcat(tx.stokesStr,'(:,1)'))), mean(eval(strcat(tx.stokesStr,'(:,2)'))), ...
        mean(eval(strcat(tx.stokesStr,'(:,3)')))];
    mlc(i1,1:3) = mlc(i1,1:3)./(norm(mlc(i1,1:3))); 
end

% maximum likelihood
for i1 = 1:length(tx.stokesNm)
    decnorm = tx.stokesNm(i1,:)*mlc';
    [vals(i1),tx.stokesInds(i1)] = max(decnorm);
end

% plot Stokes vectors on Poincare Sphere
if 0
    figure(1)
    [x,y,z]=sphere;
    surf(x,y,z,'FaceAlpha',.1,'EdgeColor','none');hold on
    plot3(normSyms(1,:),normSyms(2,:),normSyms(3,:),'.')
    axis equal
end

%% Perform cross correlation of sequences

tx.corrseq = zeros([1,length(normSyms)]);
tx.szperms = size(tx.symbPerms);
tx.dind = zeros([1,tx.szperms(1)]);
tx.dval = zeros([1,tx.szperms(1)]);
for i1 = 1:tx.szperms(1)
    for i2 = 1:tx.szperms(2)
        setInds=find(i2==tx.stokesInds);
        tx.corrseq(setInds) = tx.symbPerms(i1,i2);
    end
    % perform an cross correlation between voltage sequence and stoke sequence
    [tx.dval(i1),tx.dind(i1)] = max(xcorr(tx.corrseq,voltGroups)./(sqrt(sum(abs(voltGroups').^2)*sum(abs(tx.corrseq').^2))));
end

[tx.delayVal,tx.permInd] = max(tx.dval);
dval=tx.dval;
symDelay = tx.dind(tx.permInd) - length(voltGroups);

i1=tx.permInd;
    for i2 = 1:tx.szperms(2)
        setInds=find(i2==tx.stokesInds);
        tx.corrseq(setInds) = tx.symbPerms(i1,i2);
    end
    voltGroups_crot = tx.corrseq;
clc
end
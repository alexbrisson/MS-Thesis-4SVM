function [EVM,mlc] = SVM_evm(normSyms,M)

% This function computes the Error Vector Magnitude of Stokes
% Constellations

%% Inputs: 
% normSyms: received normalized Stokes vectors
% M: Constellation cardinality 

szst = size(normSyms);
if szst(1) == 3
normSyms = normSyms';
end

%% Seperate out groups of points and find the averages, normalize

Z = linkage(normSyms,'ward'); % linkage: Agglomerative hierarchical cluster tree for clustering
c = cluster(Z,'Maxclust',M); % clusters groups into M groups.

for i1 = 1:M
    tx.stokesStr = strcat('tx.stokes_g',num2str(i1));
    gInds = find(c==i1);
    eval([tx.stokesStr,'=normSyms(gInds,:)'])
end

for i1 = 1:M
    tx.stokesStr = strcat('tx.stokes_g',num2str(i1));
    mlc(i1,1:3) = [mean(eval(strcat(tx.stokesStr,'(:,1)'))), mean(eval(strcat(tx.stokesStr,'(:,2)'))), ...
        mean(eval(strcat(tx.stokesStr,'(:,3)')))];
    mlc(i1,1:3) = mlc(i1,1:3)./(norm(mlc(i1,1:3))); 
end

%% EVM
for i1 = 1:M
    for i2 = 1:length(eval(strcat('tx.stokes_g',num2str(i1))))
        sig_i2(i2) = norm(eval(strcat('tx.stokes_g',num2str(i1),'(i2,:)'))-mlc(i1,:))^2;
    end
    sig_i1(i1) = sum(sig_i2);
end

EVM = sqrt((1/length(normSyms))*sum(sig_i1));


end


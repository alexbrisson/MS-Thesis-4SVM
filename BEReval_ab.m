function [BER,nbits,nerrs] = BEReval_ab(nominalstokes,symbols_cr,M,txbits)
% compute decision norms to yield BER

%% Inputs:
% nominalstokes: nominal stokes vectors with M columns
% symbols_cr: counter-rotated received symbols
% M: Cardinality 
% txbits: transmitted bit sequence

bits = dec2bin(0:M-1);
rxbits = zeros([length(symbols_cr),log2(M)]);

for i = 1:length(symbols_cr)
    decNorms = symbols_cr(:,i)'*nominalstokes';
    [~,decision] = max(decNorms);
    Output=char(num2cell(bits(decision,:)));
    Output=reshape(str2num(Output),1,[]);
    rxbits(i,:) = Output;
end

rxbits = reshape(rxbits',[],1);
nbits = length(rxbits);
hammingDist = sum(rxbits ~= txbits);
nerrs=hammingDist;
BER = hammingDist/nbits;

end
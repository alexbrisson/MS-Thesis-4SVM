function txbits = transbitsGen_ab(prss,M)


bits = dec2bin(0:M-1);
txbits = zeros([length(prss),log2(M)]);

for n = 1:M
        inds=find(prss==n);
        Output=char(num2cell(bits(n,:)));
        Output=reshape(str2num(Output),1,[]);
        for i = 1:length(inds)
        txbits(inds(i),:) = Output;
        end
end

txbits = reshape(txbits',[],1);

end
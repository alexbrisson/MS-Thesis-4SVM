function stokes = jones2stokes_ab(jones)


szjones = size(jones);
if szjones(1) ~= 2
    jones = jones.';
    szjones = size(jones);
end

% pauli matrices 
sig1 = [1,0;0,-1];
sig2 = [0,1;1,0];
sig3 = [0,-1i;1i,0];

% allocate memory for stoke components
s1 = zeros([1,szjones(2)]);
s2 = zeros([1,szjones(2)]);
s3 = zeros([1,szjones(2)]);
for i=1:szjones(2)
    s1(i) = conj(jones(:,i).')*sig1*jones(:,i);
    s2(i) = conj(jones(:,i).')*sig2*jones(:,i);
    s3(i) = conj(jones(:,i).')*sig3*jones(:,i);
end

stokes = [s1;s2;s3];
function [xi,jonesnk,stokesnk,vdack] = voltOpt(stokesn,jonesnk_1,vdack_1,step)
%Steepest Ascent Algorithm - used to optimize drive voltages
%delivered to the DP-MZM for given M-ary SVM. 

% By Alex Brisson - Optical Communications Group, MSU
% 03/24

%% Inputs: 
% stokesn: Cluster averaged Stokes vectors
% jonesnk_1: k-1-th iteration Jones vectors
% vdack_1: k-1-th iteration gamma values
% step: Optimization step size

%% Ensure dimensionality of vectors is correct

szst = size(stokesn);
if szst(1) == 4
stokesn = stokesn';
szst = size(stokesn);
end

szjn = size(jonesnk_1);
if szjn(1) > 2
    jonesnk_1 = jonesnk_1.';
    szjn = size(jonesnk_1);
end
   
%% Compute Gram Matrix & current iteration cost function

for i = 1:3
    s(:,i) = stokesn(:,i)-stokesn(:,4); % col matrix
end

s = s'; % row matrix...

% Gram Matrix
Gram = s*s';
xi = det(Gram);

%% Perform Stokes 2 Jones Conversion and evaluate partial derivatives

% stokes 2 jones
jonesn = stokes2jones_ab(stokesn);

% allocate partial derivative components of G
pol1 = zeros([3,3]);
pol2 = zeros([3,3]);

% FM Index
beta = pi/2;

% Partial Derivative w.r.t Gram Matrix
for k = 1:4
    parsymb1 = strcat('dG.s',num2str(k),'1'); % xpol
    parsymb2 = strcat('dG.s',num2str(k),'2'); % ypol
for i = 1:3
    for j = 1:3
        if i == j
            pol1(i,j) = (-8*(kronDel(i,k)*conj(jonesn(:,4).')*jonesn(:,i)*jonesn(1,4)) - 8*(kronDel(4,k)*conj(jonesn(:,i).')*jonesn(:,4)));
            pol2(i,j) = (-8*(kronDel(i,k)*conj(jonesn(:,4).')*jonesn(:,i)*jonesn(2,4)) - 8*(kronDel(4,k)*conj(jonesn(:,i).')*jonesn(:,4)));
        else
            pol1(i,j) = (4*(kronDel(i,k)*conj(jonesn(:,j).')*jonesn(:,i)*jonesn(1,j)) + 4*(kronDel(j,k)*conj(jonesn(:,i).')*jonesn(:,j)*jonesn(1,i)) ...
                - 4*(kronDel(i,k)*conj(jonesn(:,4).')*jonesn(:,i)*jonesn(1,4)) - 4*(kronDel(4,k)*conj(jonesn(:,i).')*jonesn(:,4)*jonesn(1,i)) ...
                - 4*(kronDel(j,k)*conj(jonesn(:,4).')*jonesn(:,j)*jonesn(1,4)) - 4*(kronDel(4,k)*conj(jonesn(:,j).')*jonesn(:,4)*jonesn(1,j)));
            pol2(i,j) = (4*(kronDel(i,k)*conj(jonesn(:,j).')*jonesn(:,i)*jonesn(2,j)) + 4*(kronDel(j,k)*conj(jonesn(:,i).')*jonesn(:,j)*jonesn(2,i)) ...
                - 4*(kronDel(i,k)*conj(jonesn(:,4).')*jonesn(:,i)*jonesn(2,4)) - 4*(kronDel(4,k)*conj(jonesn(:,i).')*jonesn(:,4)*jonesn(2,i)) ...
                - 4*(kronDel(j,k)*conj(jonesn(:,4).')*jonesn(:,j)*jonesn(2,4)) - 4*(kronDel(4,k)*conj(jonesn(:,j).')*jonesn(:,4)*jonesn(2,j)));
        end
    end
end

    % Chain Rule 
    re_pol1=real(pol1).*beta*cos(beta*vdack_1(k,1));
    im_pol1=imag(pol1).*beta*cos(beta*vdack_1(k,2));
    re_pol2=real(pol2).*beta*cos(beta*vdack_1(k,3));
    im_pol2=imag(pol2).*beta*cos(beta*vdack_1(k,4));

    eval([parsymb1,'=re_pol1+1i.*im_pol1']);
    eval([parsymb2,'=re_pol2+1i.*im_pol2']);

end

% compute the partial xis
Ginv = inv(Gram);
parxi = zeros([3,2]);
for i = 1:4
    parsymb1 = strcat('dG.s',num2str(i),'1'); % xpol
    parsymb2 = strcat('dG.s',num2str(i),'2'); % ypol
    parxi(i,1) = xi*trace(Ginv*eval(parsymb1));
    parxi(i,2) = xi*trace(Ginv*eval(parsymb2));
end

for k = 1:4
    re_parxi1(k)=real(parxi(k,1));
    im_parxi1(k)=imag(parxi(k,1));
    re_parxi2(k)=real(parxi(k,2));
    im_parxi2(k)=imag(parxi(k,2));
end

% Recursion 
vdack(1:4,1) = vdack_1(1:4,1) + step.*(re_parxi1');
vdack(1:4,2) = vdack_1(1:4,2) + step.*(im_parxi1');
vdack(1:4,3) = vdack_1(1:4,3) + step.*(re_parxi2');
vdack(1:4,4) = vdack_1(1:4,4) + step.*(im_parxi2');

% generate kth iteration Jones Vectors
jonesk(:,1:4) = jonesn(:,1:4) + step.*parxi.';

% Normalize
for i = 1:4
jonesnk(:,i) = jonesk(:,i)./norm(jonesk(:,i));
end

% update kth iteration stokes vector
stokesnk = jones2stokes_ab(jonesnk.');

end

% Kronecker Delta Function
function c = kronDel(a,b)

if a == b
    c = 1;
else
    c = 0;
end
end


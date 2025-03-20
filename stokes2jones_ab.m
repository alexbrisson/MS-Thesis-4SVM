function [jv] = stokes2jones_ab(stokes)
% stokes 2 jones

szstokes = size(stokes);
if szstokes(1) ~= 3 % take as column matrix
    stokes = stokes';
    szstokes = size(stokes);
end


for i = 1:szstokes(2)
s1 = stokes(1,i);
s2 = stokes(2,i);
s3 = stokes(3,i);
theta = acos(s1);
cosphi = (s2/sqrt(1-(cos(theta)^2)));
sinphi = (s3/sqrt(1-(cos(theta)^2)));
if isnan(sinphi) || isinf(sinphi)
    sinphi=0;
end
if isnan(cosphi) || isinf(cosphi)
    cosphi = 0;
end
jv(:,i) = [cos(theta/2); (sqrt(1-(cos(theta/2)^2))*cosphi) + 1i*(sqrt(1-(cos(theta/2)^2))*sinphi)];
end
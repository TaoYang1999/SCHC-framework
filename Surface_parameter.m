function [sigma eta R m2] = Surface_parameter(x,z)
m0 = 0;
m2 = 0;
m4 = 0;
N = length(z);
dz = zeros(1,N);
for i=1:N-1
    s = x(i+1)-x(i);
    dz(i) = (z(i+1)-z(i))/s;
end
d2z = zeros(1,N-2);
for i=1:N-2
    s = x(i+1)-x(i);
    d2z(i) = (dz((i+1))-dz((i)))/s;
end
for i=1:N-2
    m0 = m0 + z(i)^2;
    m2 = m2 + dz(i)^2;
    m4 = m4 + d2z(i)^2; 
end
m0 = m0/length(z);
m2 = m2/length(dz);
m4 = m4/length(d2z);
% sigma = sqrt(m0)
alpha = m0*m4/m2^2;
sigma = sqrt((1-0.8968/alpha)*m0);
eta = (m4/m2)*(1/(6*pi*sqrt(3)));
R = 0.375*(pi/m4)^0.5;
end
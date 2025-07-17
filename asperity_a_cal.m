function a = asperity_a_cal(R, w)
% Input: R - radius of asperity
%        w - compression depth
% Output: a - asperity contact area

global E_roughness H kZ

E = E_roughness;

we = (pi*kZ*H/2/E)^2*R;

if 0 < w  &&  w < we

a(1) = pi*R*w;
a(2) = 1;

elseif  we < w  &&  w < 6*we

a(1) = 0.93 * pi * R^0.864 * (pi*kZ*H/2/E)^-0.272 * w^1.136;
a(2) = 2;


elseif  6*we < w  &&  w < 110*we

a(1) = 0.94 * pi * R^0.854 * (pi*kZ*H/2/E)^-0.292 * w^1.146;
a(2) = 3;


elseif w > 110*we

a(1) = 2*pi*R*w;
a(2) = 4;


end




end

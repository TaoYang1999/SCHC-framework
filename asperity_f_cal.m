function f = asperity_f_cal(R, w)
% Input: R - radius of asperity
%        w - compression depth
% Output: a - asperity contact force

global E_roughness H kZ

E = E_roughness;

we = (pi*kZ*H/2/E)^2*R;

if 0 < w  &&  w < we

f = 4/3*E*R^0.5*w^1.5;

elseif  we < w  &&  w < 6*we

f = 1.37*(pi*kZ*H/2/E)^0.15*E*R^0.575*w^1.425;

elseif  6*we < w  &&  w < 110*we

f = 1.87*(pi*kZ*H/2/E)^0.474*E*R^0.737*w^1.263;


elseif w > 110*we

f = 2*pi*R*w*H;

end


end

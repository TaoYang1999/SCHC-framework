function [K_e, F_Patch] = Constraint_conduct(side_down, F_Patch, K_Patch )



Disp = [0 0];  
side_choose_disp = side_down;
F_Patch = Force_ele_bt(side_choose_disp, Disp, F_Patch, K_Patch);
K_e = KZ_ele_bt(side_choose_disp, Disp, F_Patch, K_Patch);



end

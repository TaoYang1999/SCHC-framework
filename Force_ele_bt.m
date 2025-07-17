function Force_ele = Force_ele_bt(side_choose_disp, Disp, Force_ele, KZ)

global bt


for i=1:length(side_choose_disp)                      
    cons_DOF_choose(i) = [side_choose_disp(i)*2-1];     % x     
    cons_DOF_choose2(i) = [side_choose_disp(i)*2];      % y
end

Disp_choose = ones(length(cons_DOF_choose),1)*Disp(1);
Disp_choose2 = ones(length(cons_DOF_choose2),1)*Disp(2);

for i=1:length(cons_DOF_choose)
    if Disp(1) ~= 1    % x
        KZ(cons_DOF_choose(i),cons_DOF_choose(i)) = bt*KZ(cons_DOF_choose(i),cons_DOF_choose(i));
        Force_ele(cons_DOF_choose(i)) = Disp_choose(i)*KZ(cons_DOF_choose(i),cons_DOF_choose(i));
    end

    if Disp(2) ~= 1    % y
        KZ(cons_DOF_choose2(i),cons_DOF_choose2(i)) = bt*KZ(cons_DOF_choose2(i),cons_DOF_choose2(i));
        Force_ele(cons_DOF_choose2(i)) = Disp_choose2(i)*KZ(cons_DOF_choose2(i),cons_DOF_choose2(i));
    end
end

end
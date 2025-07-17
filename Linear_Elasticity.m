function [K, F, EA] = Linear_Elasticity(node, element, D)

global Parameter DOF

t = Parameter.t;

ele_num = length(element(:,1)); 
node_num = length(node(:,1));   

% K = zeros( 2 * node_num  , 2 * node_num  );
F = zeros( 2 * node_num  , 1 );


ele_DOF = zeros(ele_num, DOF*3);
for i=1:ele_num   
    for j=1:length(element(i,:))
        ele_DOF(i,1+2*(j-1):2+2*(j-1)) = 2*(element(i,j)-1)+1 : 2*(element(i,j)-1)+2;
    end
end

K = zeros(node_num*2, node_num*2);
for i=1:ele_num
    EA = element_area(element, node);
    ke = ele_stiff_matrix(element, node, i, EA, D, t); 
    K(ele_DOF(i,:), ele_DOF(i,:)) = K(ele_DOF(i,:), ele_DOF(i,:))+ke;
end




end
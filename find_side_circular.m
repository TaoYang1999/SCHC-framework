function side_circular = find_side_circular(element, node,X, Y, r)
ele_num = length(element(:,1));  
node_num = length(node(:,1));   
num1 = 0;
er = 1e-6;    %  *sqrt(min(element_area(element, node)));

for i=1:node_num
    
    if sqrt(abs((node(i,1)-X)^2 + (node(i,2)-Y)^2 - r^2)) <= er && node(i,2) < 2e3
        num1 = num1+1;
        side_circular(num1) = i;   
    end
   
end

end
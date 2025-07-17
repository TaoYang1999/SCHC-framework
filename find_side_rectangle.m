function [side_left, side_right, side_up, side_down] = find_side_rectangle(element, node)
ele_num = length(element(:,1)); 
node_num = length(node(:,1));   
num1 = 0;
num2 = 0;
num3 = 0;
num4 = 0;
er = 1e-2*sqrt(min(element_area(element, node)));

for i=1:node_num
    x_left = min(node(:,1));
    x_right = max(node(:,1));
    y_up = max(node(:,2));
    y_down = min(node(:,2));
    if abs(node(i,1) - x_left) <= er
        num1 = num1+1;
        side_left(num1) = i;    
    end
    if abs(node(i,1) - x_right) <= er
        num2 = num2+1;
        side_right(num2) = i;     
    end
    if abs(node(i,2) - y_up) <= er
        num3 = num3+1;
        side_up(num3) = i;    
    end
    if abs(node(i,2) - y_down) <= er
        num4 = num4+1;
        side_down(num4) = i;  
    end
end

end

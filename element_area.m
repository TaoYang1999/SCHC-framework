function A = element_area(element, node)
A = zeros(1,length(element(:,1)));
for i=1:length(element(:,1))
    x1 = node(element(i,1),1);
    x2 = node(element(i,2),1);
    x3 = node(element(i,3),1);
    y1 = node(element(i,1),2);
    y2 = node(element(i,2),2);
    y3 = node(element(i,3),2);
    A(i) = (x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2)/2;  
end
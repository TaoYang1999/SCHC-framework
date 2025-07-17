function [ ] = Plot_stress(Element,Node_s,sgm,type)


hold on
if type ~= 4
PATCH(Element, Node_s, sgm(type,:));
else
    SGM = sqrt( sgm(1,:).^2  + sgm(2,:).^2 - sgm(1,:).*sgm(2,:) + 3.*sgm(3,:).^2 );
    PATCH(Element, Node_s, SGM);
end
% patch('Faces', Element_s, 'Vertices', Node_s, 'FaceColor', 'none')
colorbar;
colormap(flipud(othercolor('RdYlBu11'))); 

axis equal
colorbar


end
function [Nnode,Nele]=MeshPlot(Node,Element)
%plot mesh
% figure
hold on;
patch('Faces',Element,'Vertices',Node,'facecolor', [.9, .9, .9], 'EdgeColor', [1, 0.5, 0.5])

axis equal
set(gcf,'color','white')
title('Mesh');
axis([-0.03 0.03 -0.025 0.01]*1)
box on


% title('Undeformed Geometry');

% num=0;
% for i=1:numel(Node(:,1))
%     num=num+1;
%     %     plot(Node(i,1),Node(i,2),'bo')
%     text(Node(i,1),Node(i,2),0,num2str(num),'Color','r');
% end
% Nnode=size(Node,1);
% Nele=size(Element,1);








function SGM = PATCH(element, node, sgm)

num_nodes = size(node, 1);

SGM = zeros(num_nodes, 1); 

node_count = zeros(num_nodes, 1); 


for ie = 1:length(element)
    nodes = element(ie, :); 
    SGM(nodes) = SGM(nodes) + sgm(ie); 
    node_count(nodes) = node_count(nodes) + 1; 
end

SGM = SGM ./ node_count;


for ie = 1:length(element)
    nodes = element(ie, :);
    x = node(nodes, 1); 
    y = node(nodes, 2);
    c = SGM(nodes); 

    patch(x, y, c, 'EdgeColor', 'interp', 'FaceColor', 'interp');
end

colorbar;
colormap(othercolor('RdYlBu11'));
set(gcf,'color','white')
axis equal;
box on



end





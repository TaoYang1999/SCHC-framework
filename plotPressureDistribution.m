function plotPressureDistribution(Cor_x, Pressure)

figure('Position', [150, 150, 1200, 500]*1); % [x, y, width, height]

num_steps = size(Pressure, 1);

y = 1:num_steps;

[X, Y] = meshgrid(Cor_x, y);

waterfall(X, Y, Pressure);



set(gca, 'FontSize', 20); 
set(gcf, 'Color', 'white'); 
grid on;
colormap(othercolor('RdYlBu11'));
colormap(flipud(othercolor('RdYlBu11'))); 
caxis([-1, 3.5]*1e5*1)

view([-20, 40]);

end
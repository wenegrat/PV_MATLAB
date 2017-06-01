figure
for i=1:7
    switch i
        case 1
            mark = 's';
             set(gca, 'ColorOrderIndex', 1);
        case 2
            mark = 's'
             set(gca, 'ColorOrderIndex', 2);
        case 3
            mark = 's'
             set(gca, 'ColorOrderIndex', 3);
        case 4
            mark = 's'
             set(gca, 'ColorOrderIndex', 4);
        case 5
            mark = 'd'
            set(gca, 'ColorOrderIndex', 5);
        case 6
            mark = 'o'
            set(gca, 'ColorOrderIndex', 5);
        case 7
            mark = 's'
            set(gca, 'ColorOrderIndex', 5);
    
    end

     hold on
     scatter(i,i, mark, 'Filled', 'MarkerEdgeColor', 'k');
     hold off
end

l = legend('$\nabla b_o = (1f)^2$','$\nabla b_o = (2f)^2$', '$\nabla b_o = (4f)^2$', '$\nabla b_o = (6f)^2$', '$Q_o = 25$ $W m^{-2}$','$Q_o = 100$ $W m^{-2}$','$Q_o = 200$ $W m^{-2}$' , 'Location', 'NorthEastOutside');
set(l, 'Interpreter' ,'latex', 'FontSize', 20)



set(gcf, 'Color', 'w');






figure
yyaxis left
plot(output.time, -output.H, 'LineWidth', 2);
set(gca, 'ylim', [-200 0])
ylabel('$\overline{H}$ (m)');
set(gca, 'FontSize', 16);
yyaxis right
plot( output.time, output.MagB.^2./f0^4, 'LineWidth', 2)
% set(gca, 'ylim', [0 7])
ylabel('$\overline{|\nabla_h b|}^2/f^4$')
set(gca, 'FontSize', 16)
xlabel('Days');

set(gca, 'xlim', [0 45]);
grid on
set(gcf, 'Color', 'w', 'Position', [368   655   728   281]);

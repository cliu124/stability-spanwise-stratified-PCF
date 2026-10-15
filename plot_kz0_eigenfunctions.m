function plot_kz0_eigenfunctions(y, v, u, p, flag, Ri_ind, kx, f0, dpdy, U_bar, d_U_bar, dd_U_bar)
% Plot and export the original complex eigenfunctions without normalization.
% Their common eigenvector amplitude and phase are retained.

y_plot = y(:);
mode_values = [v(:), u(:), p(:)];
curve_values = [real(mode_values(:,1)), imag(mode_values(:,1)), ...
    real(mode_values(:,2)), imag(mode_values(:,2)), ...
    real(mode_values(:,3)), imag(mode_values(:,3))];

Fr = 1/sqrt(flag.Ri_bulk_list(Ri_ind));
Re = flag.Re_list(1);
Pe = Re*flag.Pr_list(1);
output_folder = fullfile(pwd,sprintf('%s_Re_%g_Pe_%g',flag.mean,Re,Pe));
if ~isfolder(output_folder)
    [created,message]=mkdir(output_folder);
    if ~created
        error('Could not create results folder %s: %s',output_folder,message);
    end
end
if flag.asymptotic_low_Pr
    method = 'asymptotic_low_Pr';
else
    method = 'asymptotic';
end
output_stem = sprintf('%s_eigenfunctions_kz0_%s_Re_%g_Pe_%g_kx_%g_Fr_%g', ...
    flag.mean, method, Re, Pe, kx, Fr);

fig = figure('Color','white');
ax = axes('Parent',fig);
hold(ax,'on');
styles = {'k-o','k--s','r-^','r--v','b-d','b--+'};
labels = {'$\mathrm{Re}[\hat{v}]$','$\mathrm{Im}[\hat{v}]$', ...
    '$\mathrm{Re}[\hat{u}]$','$\mathrm{Im}[\hat{u}]$', ...
    '$\mathrm{Re}[\hat{p}]$','$\mathrm{Im}[\hat{p}]$'};
for curve_ind = 1:6
    marker_fraction = ((0:5)+(curve_ind-0.5)/6)/6;
    marker_indices = unique(1+round(marker_fraction*(numel(y_plot)-1)));
    plot(ax,curve_values(:,curve_ind),y_plot,styles{curve_ind}, ...
        'LineWidth',1.5,'MarkerSize',8,'MarkerIndices',marker_indices, ...
        'DisplayName',labels{curve_ind});
end
xlabel(ax,'Eigenfunction value','Interpreter','latex','FontSize',28);
ylabel(ax,'$y$','Interpreter','latex','FontSize',28);
legend(ax,'show','Interpreter','latex','Location','southeast', ...
    'BackgroundAlpha',0.5,'FontSize',24);
set(ax,'FontName','Times New Roman','FontSize',24,'Box','on');
ylim(ax,[-10,10]);
xlim(ax,[-0.2,0.2]);
ax.XLabel.FontSize=28;
ax.YLabel.FontSize=28;
plot_title = ['$Re=',latex_scientific(Re), ...
    ',\quad Pe=',latex_scientific(Pe), ...
    ',\quad Fr=',sprintf('%g',Fr), ...
    ',\quad k_x=',sprintf('%.4f',kx),',\quad k_z=0$'];
title(ax,plot_title,'Interpreter','latex','FontSize',24);
grid(ax,'on');
exportgraphics(fig,fullfile(output_folder,[output_stem '.png']),'Resolution',300, ...
    'Units','pixels','Padding',100,'BackgroundColor','white');

fig_f0 = figure('Color','white');
ax_f0 = axes('Parent',fig_f0);
hold(ax_f0,'on');
f0_curves = [real(f0(:)), imag(f0(:)), dd_U_bar(:)];
f0_labels = {'$\mathrm{Re}[f_0]$','$\mathrm{Im}[f_0]$','$\overline{U}^{\prime\prime}$'};
for curve_ind = 1:3
    plot(ax_f0,f0_curves(:,curve_ind),y_plot,styles{curve_ind}, ...
        'LineWidth',1.5,'MarkerSize',8, ...
        'MarkerIndices',unique(round(linspace(1,numel(y_plot),6))), ...
        'DisplayName',f0_labels{curve_ind});
end
xlabel(ax_f0,'Profile value','Interpreter','latex','FontSize',28);
ylabel(ax_f0,'$y$','Interpreter','latex','FontSize',28);
legend(ax_f0,'show','Interpreter','latex','Location','best', ...
    'BackgroundAlpha',0.5,'FontSize',24);
set(ax_f0,'FontName','Times New Roman','FontSize',24,'Box','on');
ylim(ax_f0,[-10,10]);
title(ax_f0,plot_title,'Interpreter','latex','FontSize',24);
grid(ax_f0,'on');
exportgraphics(fig_f0,fullfile(output_folder,[output_stem '_f0.png']),'Resolution',300, ...
    'Units','pixels','Padding',100,'BackgroundColor','white');

csv_table = array2table([y_plot,curve_values,real(f0(:)),imag(f0(:)), ...
    real(dpdy(:)),imag(dpdy(:)),U_bar(:),d_U_bar(:),dd_U_bar(:)],'VariableNames', ...
    {'y','v_real','v_imag','u_real','u_imag','p_real','p_imag', ...
    'f0_real','f0_imag','dpdy_real','dpdy_imag','U_bar','d_U_bar','dd_U_bar'});
writetable(csv_table,fullfile(output_folder,[output_stem '.csv']));
end

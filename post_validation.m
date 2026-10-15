function validation=post_validation(flag,output)
% Plot/export validation results and return reference data reused by asymptotics.
validation=struct();

if strcmp(flag.post,'eig_stratified_FFLWL_figure5a_validation')

%%This is the validation for spanwise stratification
 %%Comparing results of figure 5(a) of 
 %%Facchini G, Favier B, Le Gal P, Wang M, Le Bars M. The linear instability of the stratified plane Couette flow. Journal of Fluid Mechanics. 2018 Oct;853:205-34.
 %%This is for infinite Sc number, plot the growth rate varying over Re
 %%and Froude number
    data{1}.x=flag.Re_list;
    data{1}.y=output.eig_real_max_kx_kz;
    omega_FFLWL_figure5a=[665.86	0
        977.303	0.00580153
        1172.3	0.0146565
        1424.14	0.0210687
        1739.28	0.029313
        2789.12	0.0430534
        4010.41	0.0500763
        6077.53	0.0560305
        10019.5	0.060916
        20074.2	0.0650382
        50102.6	0.0677863
        ];
    data{2}.x=omega_FFLWL_figure5a(:,1);
    data{2}.y=omega_FFLWL_figure5a(:,2);
    plot_config.label_list={1,'$Re$','$Im(\omega)$'};
    plot_config.name=['validation_stratification_spanwise_eig_FFLWL_figure5a','.png'];
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
    plot_config.Markerindex=3;   
    plot_config.loglog=[0,0];
    plot_config.ylim_list=[1,-0.001,0.08];
    plot_config.xlim_list=[1,0,50200];
    plot_config.print_size=[1,1000,900];
    plot_config.xtick_list=[1,10000,20000,30000,40000,50000];
    plot_config.ytick_list=[1,0,0.02,0.04,0.06,0.08];
    plot_config.user_color_style_marker_list={'b*','ko'};
    plot_line(data,plot_config);
    

elseif strcmp(flag.post,'eig_stratified_FFLWL_figure10_validation')

%%This is the validation for spanwise stratification
     %%Comparing results of figure 10(b) of 
     %%Facchini G, Favier B, Le Gal P, Wang M, Le Bars M. The linear instability of the stratified plane Couette flow. Journal of Fluid Mechanics. 2018 Oct;853:205-34.
     %%This is for a finite Sc number, and only for one parameter
    data{1}.x=real(output.eig_val{1}/(-1i));
    data{1}.y=imag(output.eig_val{1}/(-1i));
    plot_config.label_list={1,'$Re(\omega)$','$Im(\omega)$'};
    plot_config.xlim_list=[1,-2,2];
    plot_config.ylim_list=[1,-1,0.05];
    plot_config.print_size=[1,1100,900];
    plot_config.loglog=[0,0];
    plot_config.xtick_list=[1,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2];
    plot_config.ytick_list=[1,-1,-0.8,-0.6,-0.4,-0.2,0];
    plot_config.name=['validation_stratification_spanwise_eig_FFLWL_figure10a','.png'];
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
    plot_config.Markerindex=3;
    plot_config.user_color_style_marker_list={'ksquare','rsquare'};
    plot_line(data,plot_config);
    
    
    data{1}.x=real(output.eig_val{1}/(-1i));
    data{1}.y=imag(output.eig_val{1}/(-1i));
    omega_FFLWL_figure10b=[0	0.00793835
        -0.336207	-0.10153
        0.318966	-0.102134
        -0.00862069	-0.122127
        -0.491379	-0.167808
        0.482759	-0.168706];
    data{2}.x=omega_FFLWL_figure10b(:,1);
    data{2}.y=omega_FFLWL_figure10b(:,2);
    plot_config.xlim_list=[1,-0.6,0.6];
    plot_config.ylim_list=[1,-0.2,0.075];
    plot_config.print_size=[1,600,900];
    plot_config.xtick_list=[1,-0.5,0,0.5];
    plot_config.ytick_list=[1,-0.2,-0.15,-0.1,-0.05,0,0.05];
    plot_config.name=['validation_stratification_spanwise_eig_FFLWL_figure10b','.png'];
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
    plot_config.user_color_style_marker_list={'k*','rsquare'};
    plot_line(data,plot_config);

elseif strcmp(flag.post,'eig_stratified_spanwise_LHBLCF_figure6b_validation')
    data{1}.x=flag.kx_list;
    data{1}.y=flag.kz_list;
    data{1}.z=output.eig_real_max{1}';
    data{1}.z(find(data{1}.z<0))=NaN;

    data{2}.x=1*ones(20,1);
    data{2}.y=linspace(0,19,20);
    plot_config.colormap='parula';
    plot_config.label_list={1,'$k_x$','$k_z$'};
    plot_config.xlim_list=[1,0,2.6];
    plot_config.ylim_list=[1,0,19];
    plot_config.zlim_list=[1,0,0.05];
    plot_config.ztick_list=[1,0.01,0.02,0.03,0.04,0.05];
    plot_config.print_size=[1,1100,900];
    plot_config.loglog=[0,0];
    plot_config.xtick_list=[1,0,0.5,1,1.5,2,2.5];
    plot_config.ytick_list=[1,0,2,4,6,8,10,12,14,16,18];
    plot_config.name=['validation_stratification_spanwise_eig_LHBLCF_figure6b','.png'];
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
    plot_config.user_color_style_marker_list={'k--'};
    plot_contour(data,plot_config);
    
    
elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_DCB_figure4a') || ...
        (strcmp(flag.post,'eig_stratified_spanwise_tanh_asymptotic') && strcmp(flag.mean,'tanh'))
    validation.DCB_figure4a_Fh1=[0	0.189744
            0.238569	0.178917
            0.357853	0.164103
            0.437376	0.150427
            0.516899	0.135613
            0.55666	0.120798
            0.636183	0.107123
            0.675944	0.0917379
            0.715706	0.0774929
            0.755467	0.0638177
            0.795229	0.048433
            0.795229	0.0353276
            0.83499	0.019943
            0.83499	0.00740741
            ];
    validation.DCB_figure4a_Fh01=[0	0.189714
            0.874751	0.188
            1.5507	0.184571
            2.18688	0.18
            2.74354	0.174286
            3.26044	0.168
            3.77734	0.160571
            4.334	0.151429
            4.77137	0.143429
            5.36779	0.131429
            5.88469	0.12
            6.36183	0.108
            6.64016	0.101143
            7.27634	0.0817143
            7.71372	0.0657143
            7.99205	0.0537143
            8.19085	0.0451429
            8.54871	0.0268571
            8.66799	0.0188571
            ];
    
    validation.DCB_figure4a_Fh005=[0	0.189714
            1.0338	0.189143
            2.06759	0.187429
            3.06163	0.185143
            4.09543	0.181143
            5.08946	0.176
            6.12326	0.170857
            7.1173	0.163429
            8.19085	0.156
            9.18489	0.146857
            10.1789	0.137143
            11.2127	0.126286
            12.2465	0.114286
            13.3201	0.1
            14.2744	0.0851429
            15.3082	0.0685714
            16.3419	0.0468571
            17.336	0.0182857
            ];
    
    data{1}.x=flag.kz_list;
    data{1}.y=output.eig_real_max{1};
    data{2}.x=flag.kz_list;
    data{2}.y=output.eig_real_max{2};
    data{3}.x=flag.kz_list;
    data{3}.y=output.eig_real_max{3};
    if strcmp(flag.post,'eig_stratified_spanwise_tanh_asymptotic')
        for plot_Ri_ind=1:3
            data{plot_Ri_ind}.x=flag.Fh_kz_list*sqrt(flag.Ri_bulk_list(plot_Ri_ind));
        end
    end
    data{4}.x=validation.DCB_figure4a_Fh1(:,1);
    data{4}.y=validation.DCB_figure4a_Fh1(:,2);
    data{5}.x=validation.DCB_figure4a_Fh01(:,1);
    data{5}.y=validation.DCB_figure4a_Fh01(:,2);
    data{6}.x=validation.DCB_figure4a_Fh005(:,1);
    data{6}.y=validation.DCB_figure4a_Fh005(:,2);
    % Apply the same tanh Fr=1 cutoff to numerical and reference validation curves.
    for cutoff_curve_ind=1:numel(data)
        cutoff_Fr=1/sqrt(flag.Ri_bulk_list(mod(cutoff_curve_ind-1,3)+1));
        if abs(cutoff_Fr-1)<1e-12
            cutoff_keep=data{cutoff_curve_ind}.x<=0.82;
            data{cutoff_curve_ind}.x=data{cutoff_curve_ind}.x(cutoff_keep);
            data{cutoff_curve_ind}.y=data{cutoff_curve_ind}.y(cutoff_keep);
        end
    end
    plot_config.Markerindex=3;
    plot_config.loglog=[0,0];
    plot_config.ylim_list=[1,0.02,0.2];
    plot_config.ytick_list=[1,0,0.1,0.2];
    plot_config.xlim_list=[1,0,20];
    plot_config.xtick_list=[1,0,5,10,15,20];
    plot_config.user_color_style_marker_list={'k-','r--','b-.','k^','rsquare','bo'};
    plot_config.label_list={1,'$k_z$','$\sigma$'};
    plot_config.name=['eig_stratified_spanwise_validation_stratification_z_tanh_DCB_figure4a','.png'];
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
    plot_line(data,plot_config);

elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_DCB_figure5')
    DCB_figure5_Fh1=[0	0.143026
        0.159057	0.142711
        0.324006	0.141139
        0.488954	0.138939
        0.653903	0.135481
        0.818851	0.131081
        0.9838	0.126051
        1.14875	0.120079
        1.3137	0.113477
        1.46686	0.106248
        1.63181	0.098389
        1.80265	0.0902161
        2.02062	0.0814145
        2.10898	0.0782711
        2.24448	0.074499
        2.4271	0.0700982
        2.56848	0.0669548
        2.72754	0.0641257
        2.95729	0.060668
        3.12813	0.0584676
        3.29897	0.0562672
        3.44624	0.0546955
        3.59941	0.0531238
        3.79971	0.0512377
        4	0.049666
        ];
    DCB_figure5_Fh01=[0	0.143026
        0.159057	0.142397
        0.324006	0.141139
        0.488954	0.138939
        0.653903	0.135796
        0.818851	0.131395
        0.9838	0.12668
        1.14286	0.121022
        1.30781	0.114735
        1.47275	0.107505
        1.6377	0.0996464
        1.79676	0.0911591
        1.96171	0.0829862
        2.12666	0.0757564
        2.29161	0.0700982
        2.45066	0.0660118
        2.61561	0.062554
        2.78056	0.0594106
        2.94551	0.0568959
        3.11046	0.0543811
        3.26951	0.0528094
        3.43446	0.0509234
        3.59352	0.0490373
        3.75847	0.048723
        3.92342	0.0455796
        ];
    
    %%x axis is the kz times Fh, that is 1/sqrt(Ri_b) in this code
    %%implementation..
    data{1}.x=flag.Fh_kz_list;%flag.kz_list/sqrt(flag.Ri_bulk_list(1));
    data{1}.y=output.eig_real_max{1};
    data{2}.x=flag.Fh_kz_list;%flag.kz_list/sqrt(flag.Ri_bulk_list(2));
    data{2}.y=output.eig_real_max{2};
    data{3}.x=flag.Fh_kz_list;%flag.kz_list/sqrt(flag.Ri_bulk_list(3));
    data{3}.y=output.eig_real_max{3};
    data{4}.x=DCB_figure5_Fh1(:,1);
    data{4}.y=DCB_figure5_Fh1(:,2);
    data{5}.x=DCB_figure5_Fh01(:,1);
    data{5}.y=DCB_figure5_Fh01(:,2);
    
    plot_config.Markerindex=3;
    plot_config.loglog=[0,0];
    plot_config.ylim_list=[1,0,0.16];
    plot_config.ytick_list=[1,0,0.04,0.08,0.12,0.16];
    plot_config.xlim_list=[1,0,4];
    plot_config.xtick_list=[1,0,1,2,3,4];
    plot_config.user_color_style_marker_list={'k-','r--','b-.','k^','rsquare','bo'};
    plot_config.label_list={1,'$Fr\,k_z$','$\sigma$'};
    plot_config.name=['eig_stratified_spanwise_validation_stratification_z_tanh_bickley_DCB_figure5','.png'];
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
    plot_line(data,plot_config);
elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_growth_rate_kx') || strcmp(flag.post,'eig_stratified_spanwise_bickley_growth_rate_kx')
    for Ri_ind=1:length(flag.Ri_bulk_list)
        for kz_ind=1:length(flag.kz_list)
            output.eig_real_max{Ri_ind};
            data{kz_ind}.x=flag.kx_list;
            data{kz_ind}.y=output.eig_real_max{Ri_ind}(:,kz_ind);
        end
        Ri=flag.Ri_bulk_list(Ri_ind);
        Fr=sqrt(1/Ri);
        plot_config.Markerindex=3;
        % plot_config.loglog=[0,0];
        plot_config.ylim_list=[1,0,0.25];
        plot_config.ytick_list=[1,0,0.05,0.1,0.15,0.2,0.25,0.3];
        if strcmp(flag.post,'eig_stratified_spanwise_tanh_growth_rate_kx')
            plot_config.xlim_list=[1,0,1.05];
            plot_config.xtick_list=[1,0,0.2,0.4,0.6,0.8,1,1.2];
        elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_growth_rate_kx')
            plot_config.xlim_list=[1,0,2];
            plot_config.xtick_list=[1,0,0.5,1,1.5,2];
        end
        plot_config.user_color_style_marker_list={'k-','r--','b-.','m:','c:'};
        plot_config.linewidth=2.5;
        plot_config.legend_list={1,'$k_z=0$','$k_z=2$','$k_z=4$','$k_z=6$','$k_z=8$'};
        % plot_config.user_color_style_marker_list={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30'};
        plot_config.label_list={1,'$k_x$','$\mathrm{Re}[\lambda]$'};
        %plot_config.name=['eig_stratified_spanwise_stratification_z_',flag.mean,'_growth_rate_kx_Fr=',num2str(Fr),'.png'];
        plot_config.name=sprintf('%s_numerical_eig_real_max_vs_kx_Re_%g_Pe_%g_Fr_%g.png', ...
            flag.mean,flag.Re_list(1),flag.Re_list(1)*flag.Pr_list(1),Fr);
        plot_config.name=fullfile(flag.results_folder,plot_config.name);

        comparison_title = ['$Re=',latex_scientific(flag.Re_list(1)), ...
            ',\quad Pe=',latex_scientific(flag.Re_list(1)*flag.Pr_list(1)), ...
            ',\quad Fr=',num2str(Fr),'$'];
        plot_config.title_list={1,comparison_title};
        % Export each plotted curve as adjacent x/y columns with headers.
        csv_curve_count = length(flag.kz_list);
        csv_export = cell(numel(flag.kx_list) + 1, 2*csv_curve_count);
        for csv_curve_ind = 1:csv_curve_count
            csv_label = sprintf('Numerical Fr=%g, kz=%g', Fr, flag.kz_list(csv_curve_ind));
            csv_col = 2*csv_curve_ind - 1;
            csv_export(1, csv_col:csv_col+1) = {[csv_label ', kx'], [csv_label ', real_lambda']};
            csv_export(2:end, csv_col) = num2cell(data{csv_curve_ind}.x(:));
            csv_export(2:end, csv_col+1) = num2cell(data{csv_curve_ind}.y(:));
        end
        [csv_folder, csv_stem] = fileparts(plot_config.name);
        writecell(csv_export, fullfile(csv_folder, [csv_stem '.csv']));

        plot_line(data,plot_config);
    end

elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_CLDM_figure1')

    % Transpose so rows correspond to the vertical axis (kz)
    growth_rate = real(output.eig_real_max{length(flag.Ri_bulk_list)}).';

    figure('Color', 'w', 'Position', [100, 100, 800, 700]);
    hold on;

    % Digitized reference: solid black lines
    for lev = 0.07:0.02:0.15
        filename = fullfile(sprintf('contour_%.2f.csv', lev));
        T = readtable(filename);

        plot(T.kx, T.kz, 'k-', 'LineWidth', 1.0, ...
            'HandleVisibility', 'off');
    end

    % Numerical results: dashed red contours
    levels = 0.01:0.01:0.15;

    [C, h] = contour(flag.kx_list, flag.kz_list, growth_rate, levels, ...
        'LineColor', [0.85, 0.1, 0.1], ...
        'LineStyle', '--', ...
        'LineWidth', 1.0, ...
        'HandleVisibility', 'off');

    clabel(C, h, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 13, ...
        'Color', [0.85, 0.1, 0.1], ...
        'LabelSpacing', 350);

    % Legend entries
    h_ref = plot(nan, nan, 'k-', 'LineWidth', 1.0);
    h_num = plot(nan, nan, '--', ...
        'Color', [0.85, 0.1, 0.1], 'LineWidth', 1.0);

    legend([h_ref, h_num], ...
        {'Digitized reference (0.07–0.15)', 'Numerical results'}, ...
        'Location', 'northeast', 'Box', 'off');

    xlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);

    % Your kz corresponds to ky in the reference figure
    ylabel('$k_y$', 'Interpreter', 'latex', ...
        'FontSize', 20, 'Rotation', 0);

    set(gca, ...
        'FontName', 'Times New Roman', ...
        'FontSize', 15, ...
        'LineWidth', 0.8, ...
        'TickDir', 'in', ...
        'Box', 'on', ...
        'XLim', [0, 2], ...
        'YLim', [0, 10], ...
        'XTick', 0:0.5:2, ...
        'YTick', 0:1:10);

    grid off;

    exportgraphics(gcf, fullfile(flag.results_folder,'Bickley_growth_rate_contours.pdf'), ...
        'ContentType', 'vector');


elseif any(strcmp(flag.post, {'eig_stratified_spanwise_tanh_asymptotic', ...
        'eig_stratified_spanwise_bickley_asymptotic'}))
    % Plot and save these modes in the asymptotic blocks below.
else
    error('Wrong flag.post')
end
end

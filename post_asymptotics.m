function post_asymptotics(flag,output)
% Plot/export asymptotic results independently of validation postprocessing.
% Use the physical grid returned by the solver, as in the main script.
flag.kz_list=output.kz_list;
scaled_kz=any(strcmp(flag.post,{'eig_stratified_spanwise_tanh_asymptotic', ...
    'eig_stratified_spanwise_bickley_asymptotic'}));

% Local digitized reference data for tanh comparisons.
validation=struct();
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


    validation.DCB_figure4b_Fh1=[0	0.19019
        0.05	0.18972
        0.1017	0.18785
        0.1517	0.18505
        0.2034	0.18084
        0.2534	0.17523
        0.2977	0.16893
        0.353	0.1616
        0.3861	0.15654
        0.405	0.15321
        0.431	0.14856
        0.4521	0.14472
        0.4741	0.13994
        0.4969	0.1355
        0.5163	0.13061
        0.5361	0.12593
        0.5564	0.12083
        0.5759	0.11589
        0.5948	0.11075
        0.6131	0.10547
        0.6309	0.10008
        0.6489	0.09447
        0.6658	0.08883
        0.6827	0.0829
        0.6976	0.07718
        0.7132	0.071
        0.7269	0.06486
        0.7407	0.05849
        0.7532	0.05195
        0.7653	0.04549
        0.7756	0.03876
        0.7857	0.03195
        0.7945	0.02519
        0.8025	0.01817
        0.8105	0.01114
        0.8138	0.00888
        ];

    validation.DCB_figure4b_Fh01=[0	0.19019
        0.0241	0.19019
        0.0759	0.18925
        0.1259	0.18692
        0.1776	0.18364
        0.2293	0.17897
        0.2793	0.17383
        0.331	0.1676
        0.381	0.16028
        0.4328	0.15187
        0.4845	0.14206
        0.5362	0.13146
        0.5862	0.12056
        0.6379	0.10794
        0.6879	0.09377
        0.7397	0.0771
        0.7914	0.05717
        0.8414	0.03302
        0.8672	0.01916
        ];

    validation.DCB_figure4b_Fh005=[0	0.19042
        0.05	0.18995
        0.1017	0.18808
        0.1517	0.18551
        0.2034	0.18154
        0.2534	0.17664
        0.3052	0.17103
        0.3569	0.16402
        0.4069	0.15607
        0.4586	0.1472
        0.5103	0.13692
        0.5603	0.1264
        0.6121	0.11449
        0.6638	0.10117
        0.7138	0.08575
        0.7655	0.06822
        0.8155	0.04696
        0.8672	0.01893
        ];


if flag.asymptotic
    for Ri_ind=1:length(flag.Ri_bulk_list)
        if scaled_kz
            % Evaluate asymptotics on the same physical kz grid as the eigensolver.
            flag.kz_list=flag.Fh_kz_list*sqrt(flag.Ri_bulk_list(Ri_ind));
        end

        %Get the eigenvector associated with kz=0. Here assume the first
        %value in flag.kz_list is 0. Thus just take index=1 for the kz.
        eig_vec_kz0=output.eig_vec{1,1,Ri_ind,1,1,1};
        eig_val_kz0=output.eig_val{1,1,Ri_ind,1,1,1};
        [val,ind]=max(real(eig_val_kz0)); %get the index with the maximal real part
        
        lambda_kz0=eig_val_kz0(ind); %get the lambda for most unstable mode
        c_kz0=lambda_kz0/(-1i*flag.kx_list); %get corresponding c valie, lambda=-1i*kx*c

        v_kz0=eig_vec_kz0(1:flag.Ny_full-2); %get the v velocity mode.
        u_kz0=-output.D1_bc*v_kz0/(1i*flag.kx_list); %solve u velocity using continuity equation with kz=0

        dpdy_kz0=-(1i*flag.kx_list*(output.U_bar-c_kz0).*v_kz0-1/output.Re*(output.D2_bc-output.kx^2*output.I_bc)*v_kz0);
        p_kz0=(1i*flag.kx_list*(output.U_bar-c_kz0).*u_kz0+v_kz0.*output.d_U_bar-1/output.Re*(output.D2_bc-output.kx^2*output.I_bc)*u_kz0)/(-1i*flag.kx_list);
        
        f0_kz0=1i*dpdy_kz0/flag.kx_list./((output.U_bar-c_kz0).^2);
        
        for y_ind=1:flag.Ny_full-2
            %to check directly performing the integration also match p_kz0
            %computed from x-momentum equation. p_kz0_int is not used. 
            p_kz0_int(y_ind)=sum(dpdy_kz0(y_ind:end).*output.y_diff_unbounded_bc.weight(y_ind:end));
        end

        int_weight=output.y_diff_unbounded_bc.weight; %integration weight.
        den=sum((dpdy_kz0.^2+flag.kx_list^2*p_kz0.^2)./(output.U_bar-c_kz0).^3.*int_weight);
        c1_const=flag.kx_list^2/2*sum(p_kz0.^2.*int_weight)/den; %formula of c1 without kz^2
        
        c1_const_list(Ri_ind)=c1_const; 

        c1_asymptotic=c1_const*flag.kz_list.^2; %c1 as a function of kz
        Fr2=1/flag.Ri_bulk_list(Ri_ind);
        c_asymptotic=c_kz0+Fr2*c1_asymptotic;
        lambda_asymptotic{Ri_ind}=-1i*flag.kx_list*c_asymptotic;

        % Plot/export raw kz=0 eigenfunctions; retain relative magnitudes.
        plot_kz0_eigenfunctions(output.y_bc,v_kz0,u_kz0,p_kz0,flag,Ri_ind,flag.kx_list(1), ...
            f0_kz0,dpdy_kz0,output.U_bar,output.d_U_bar,output.dd_U_bar);
    end

    if any(strcmp(flag.post, {'eig_stratified_spanwise_tanh_DCB_figure4a', ...
            'eig_stratified_spanwise_tanh_asymptotic', 'eig_stratified_spanwise_bickley_asymptotic'}))
    
        % Reuse the computed growth rates for each requested horizontal axis.
        plot_axis_names={'kz'};
        if scaled_kz
            plot_axis_names={'kz','Fr_kz'};
        end
        for plot_axis_ind=1:numel(plot_axis_names)
        for plot_component_ind=1:2 % Real and imaginary parts of the same eigenvalues.
                csv_value_name='real_lambda';
                if plot_component_ind==2
                    csv_value_name='imag_lambda';
                end
                csv_axis_name=plot_axis_names{plot_axis_ind};
                plot_kz_list=flag.kz_list;
                data = cell(1,6); % Discard curves from earlier plots.
                data{1}.x=plot_kz_list;
                data{1}.y=output.eig_real_max{1};
                data{2}.x=plot_kz_list;
                data{2}.y=output.eig_real_max{2};
                data{3}.x=plot_kz_list;
                data{3}.y=output.eig_real_max{3};
                % Asymptotic curves precede reference markers in plots and CSV files.
                data{4}.x=plot_kz_list;
                data{4}.y=real(lambda_asymptotic{1});
                data{5}.x=plot_kz_list;
                data{5}.y=real(lambda_asymptotic{2});
                data{6}.x=plot_kz_list;
                data{6}.y=real(lambda_asymptotic{3});
                if plot_component_ind==2
                    for plot_Ri_ind=1:3
                        data{plot_Ri_ind}.y=output.eig_imag_max{1,1,plot_Ri_ind,1};
                        data{plot_Ri_ind+3}.y=imag(lambda_asymptotic{plot_Ri_ind});
                    end
                end

                if strcmp(flag.mean,'tanh') && plot_component_ind==1 && strcmp(csv_axis_name,'kz')
                    data{7}.x=validation.DCB_figure4a_Fh1(:,1);
                    data{7}.y=validation.DCB_figure4a_Fh1(:,2);
                    data{8}.x=validation.DCB_figure4a_Fh01(:,1);
                    data{8}.y=validation.DCB_figure4a_Fh01(:,2);
                    data{9}.x=validation.DCB_figure4a_Fh005(:,1);
                    data{9}.y=validation.DCB_figure4a_Fh005(:,2);
                    csv_curve_types = {'Numerical', 'Asymptotic', 'DCB07'};
                    plot_config.user_color_style_marker_list={'k-','r-','b-','k--','r--','b--','k^','rsquare','bo'};
                elseif strcmp(flag.mean,'tanh') && plot_component_ind==1 && strcmp(csv_axis_name,'Fr_kz')
                    data{7}.x=validation.DCB_figure4b_Fh1(:,1);
                    data{7}.y=validation.DCB_figure4b_Fh1(:,2);
                    data{8}.x=validation.DCB_figure4b_Fh01(:,1);
                    data{8}.y=validation.DCB_figure4b_Fh01(:,2);
                    data{9}.x=validation.DCB_figure4b_Fh005(:,1);
                    data{9}.y=validation.DCB_figure4b_Fh005(:,2);
                    csv_curve_types = {'Numerical', 'Asymptotic', 'DCB07'};
                    plot_config.user_color_style_marker_list={'k-','r-','b-','k--','r--','b--','k^','rsquare','bo'};
                else
                    csv_curve_types = {'Numerical', 'Asymptotic'};
                    plot_config.user_color_style_marker_list={'k-','r--','b:','k-','r--','b:'};
                end
                plot_config.legend_list = {1};
                for plot_curve_ind = 1:numel(data)
                    plot_Fr = 1/sqrt(flag.Ri_bulk_list(mod(plot_curve_ind-1,3)+1));
                    plot_config.legend_list{plot_curve_ind+1} = sprintf('%s, Fr=%g', ...
                        csv_curve_types{ceil(plot_curve_ind/3)}, plot_Fr);
                end

                if scaled_kz
                    for plot_curve_ind=1:6
                        plot_Fr=1/sqrt(flag.Ri_bulk_list(mod(plot_curve_ind-1,3)+1));
                        if strcmp(csv_axis_name,'Fr_kz')
                            data{plot_curve_ind}.x=flag.Fh_kz_list;
                        else
                            % Each Fr has its own physical kz grid, shared by both methods.
                            data{plot_curve_ind}.x=flag.Fh_kz_list/plot_Fr;
                        end
                    end
                end
                % Digitized tanh reference abscissae are physical kz values.
                % if scaled_kz && strcmp(csv_axis_name,'Fr_kz')
                %     for plot_curve_ind=7:numel(data)
                %         plot_Fr=1/sqrt(flag.Ri_bulk_list(mod(plot_curve_ind-1,3)+1));
                %         data{plot_curve_ind}.x=data{plot_curve_ind}.x*plot_Fr;
                %     end
                % end
                % Discard tanh Fr=1 points beyond Fr*kz=0.82 before plots/CSV export.
                if strcmp(flag.mean,'tanh')
                    for cutoff_curve_ind=1:numel(data)
                        cutoff_Fr=1/sqrt(flag.Ri_bulk_list(mod(cutoff_curve_ind-1,3)+1));
                        if abs(cutoff_Fr-1)<1e-12
                            % For Fr=1, kz and Fr*kz are the same coordinate.
                            cutoff_keep=data{cutoff_curve_ind}.x<=0.82;
                            data{cutoff_curve_ind}.x=data{cutoff_curve_ind}.x(cutoff_keep);
                            data{cutoff_curve_ind}.y=data{cutoff_curve_ind}.y(cutoff_keep);
                        end
                    end
                end
                plot_config.Markerindex=3;
                plot_config.loglog=[0,0];
                plot_config.ylim_list=[1,0,0.2];
                plot_config.ytick_list=[1,0,0.1,0.2];
                plot_config.xlim_list=[1,0,25];
                plot_config.xtick_list=[1,0,5,10,15,20,25];
                plot_config.linewidth=2.5;
                plot_config.fontsize_legend=32;
                plot_config.label_list={1,'$k_z$','$\mathrm{Re}[\lambda]$'};
                if scaled_kz
                    if strcmp(csv_axis_name,'Fr_kz')
                        plot_config.label_list={1,'$Fr\,k_z$','$\mathrm{Re}[\lambda]$'};
                        if strcmp(flag.mean,'tanh')
                            plot_config.xlim_list=[1,0,max(flag.Fh_kz_list)];
                        else
                            plot_config.xlim_list=[1,0,4];
                        end
                        plot_config.xtick_list=[1,linspace(0,plot_config.xlim_list(3),5)];
                    else
                        plot_kz_max=max(flag.Fh_kz_list)*sqrt(max(flag.Ri_bulk_list));
                        if strcmp(flag.mean,'tanh')
                            % Leave room to the right of the curves for the legend.
                            plot_kz_max=5*ceil(1.35*plot_kz_max/5);
                        end
                        plot_config.xlim_list=[1,0,plot_kz_max];
                        plot_config.xtick_list=[1,linspace(0,plot_kz_max,5)];
                    end
                end
                plot_config.name= sprintf([flag.mean,'_growth_rate_vs_',csv_axis_name,'_asymptotic_Re_%g_Pe_%g_kx_%g.png'], ...
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), flag.kx_list(1));
                plot_config.name=fullfile(flag.results_folder,plot_config.name);
                % if scaled_kz
                %     % Distinct symbols/styles distinguish even exactly coincident curves.
                %     plot_config.user_color_style_marker_list = ...
                %         {'k-','r--','b:','k-','r--','b:'};
                % end
                if plot_component_ind==2
                    plot_config.label_list{3}='$\mathrm{Im}[\lambda]$';
                    % Imaginary parts may be negative or zero: choose their own limits.
                    plot_imag_values=cellfun(@(curve) curve.y(:),data,'UniformOutput',false);
                    plot_imag_values=vertcat(plot_imag_values{:});
                    plot_imag_values=plot_imag_values(isfinite(plot_imag_values));
                    plot_imag_min=min(plot_imag_values);
                    plot_imag_max=max(plot_imag_values);
                    plot_imag_padding=0.05*max([plot_imag_max-plot_imag_min,abs(plot_imag_min),abs(plot_imag_max),1e-6]);
                    % Reserve additional headroom for the upper-right legend.
                    plot_imag_upper_padding=max(plot_imag_padding,0.4*(plot_imag_max-plot_imag_min));
                    plot_config.ylim_list=[1,plot_imag_min-plot_imag_padding,plot_imag_max+plot_imag_upper_padding];
                    plot_config.ytick_list=[1,linspace(plot_imag_min-plot_imag_padding,plot_imag_max+plot_imag_upper_padding,5)];
                    plot_config.name=strrep(plot_config.name,'.png','_imag.png');
                else
                    plot_config.name=strrep(plot_config.name,'.png','_real.png');
                end
                plot_line(data,plot_config);
                if strcmp(flag.mean,'tanh')
                    comparison_legend=findobj(gcf,'Type','legend');
                    if strcmp(csv_axis_name,'Fr_kz')
                        set(comparison_legend,'BackgroundAlpha',0.5,'Location','southwest');
                    else
                        set(comparison_legend,'BackgroundAlpha',0.5,'Location','northeast');
                    end
                end
                comparison_title = ['$Re=',latex_scientific(flag.Re_list(1)), ...
                    ',\quad Pe=',latex_scientific(flag.Re_list(1)*flag.Pr_list(1)), ...
                    ',\quad k_x=',sprintf('%.4f',flag.kx_list(1)),'$'];
                title(gca,comparison_title,'Interpreter','latex','FontSize',40);
        
                if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
                    % Exactly three markers per curve, with six interleaved phases.
                    plot_marker_count = 10;
                    plot_markers = {'o','s','^','+','x','d'};
                    plot_handles = flipud(findobj(gca,'Type','line'));
                    for plot_curve_ind = 1:6
                        % Match plotted data rather than relying on pre-existing markers.
                        plot_match = find(arrayfun(@(h) ...
                            isequal(h.XData(:),data{plot_curve_ind}.x(:)) && ...
                            isequal(h.YData(:),data{plot_curve_ind}.y(:)),plot_handles),1);
                        if isempty(plot_match)
                            continue;
                        end
                        plot_handle = plot_handles(plot_match);
                        plot_handles(plot_match) = []; % Each curve is styled once.
                        plot_point_count = numel(plot_handle.XData);
                        plot_marker_phase = (plot_curve_ind-0.5)/6;
                        plot_marker_fraction = ((0:plot_marker_count-1)+plot_marker_phase)/plot_marker_count;
                        plot_marker_indices = 1 + round(plot_marker_fraction*(plot_point_count-1));
                        plot_handle.Marker = plot_markers{plot_curve_ind};
                        plot_handle.MarkerIndices = unique(plot_marker_indices);
                        plot_handle.MarkerSize = 10;
                        plot_handle.MarkerFaceColor = 'none';
                    end
                    % Refresh the saved figure after updating markers and legend symbols.
                end

                % Save the final title and marker settings for either base flow.
                exportgraphics(gcf,plot_config.name,'Resolution',300, ...
                    'Units','pixels','Padding',100,'BackgroundColor','white');

                % Save every plotted curve as paired kz/lambda columns for replotting.
                % Reference curves may have different lengths; unused cells stay blank.
                csv_curve_count = numel(data);
                csv_row_count = max(cellfun(@(curve) numel(curve.x), data(1:csv_curve_count)));
                csv_export = cell(csv_row_count + 1, 2*csv_curve_count);
                for csv_curve_ind = 1:csv_curve_count
                    csv_Fr = 1/sqrt(flag.Ri_bulk_list(mod(csv_curve_ind-1, 3) + 1));
                    csv_label = sprintf('%s Fr=%g', csv_curve_types{ceil(csv_curve_ind/3)}, csv_Fr);
                    csv_col = 2*csv_curve_ind - 1;
                    csv_export(1, csv_col:csv_col+1) = {[csv_label ', ' csv_axis_name], [csv_label ', ' csv_value_name]};
                    csv_export(2:numel(data{csv_curve_ind}.x)+1, csv_col) = num2cell(data{csv_curve_ind}.x(:));
                    csv_export(2:numel(data{csv_curve_ind}.y)+1, csv_col+1) = num2cell(data{csv_curve_ind}.y(:));
                end
                csv_filename = sprintf([flag.mean,'_growth_rate_vs_',csv_axis_name,'_asymptotic_Re_%g_Pe_%g_kx_%g.csv'], ...
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), flag.kx_list(1));
                if plot_component_ind==2
                    csv_filename=strrep(csv_filename,'.csv','_imag.csv');
                else
                    csv_filename=strrep(csv_filename,'.csv','_real.csv');
                end
                csv_filename=fullfile(flag.results_folder,csv_filename);
                writecell(csv_export, csv_filename);
        end % Real/imaginary comparison versions.
        end % Horizontal-axis versions; no eigenvalues or asymptotics are recomputed.
        
    end 

end


if flag.asymptotic_low_Pr
    for Ri_ind=1:length(flag.Ri_bulk_list)
        if scaled_kz
            % Evaluate asymptotics on the same physical kz grid as the eigensolver.
            flag.kz_list=flag.Fh_kz_list*sqrt(flag.Ri_bulk_list(Ri_ind));
        end

        %Get the eigenvector associated with kz=0. Here assume the first
        %value in flag.kz_list is 0. Thus just take index=1 for the kz.
        eig_vec_kz0=output.eig_vec{1,1,Ri_ind,1,1,1};
        eig_val_kz0=output.eig_val{1,1,Ri_ind,1,1,1};
        [val,ind]=max(real(eig_val_kz0)); %get the index with the maximal real part
        
        lambda_kz0=eig_val_kz0(ind); %get the lambda for most unstable mode
        c_kz0=lambda_kz0/(-1i*flag.kx_list); %get corresponding c valie, lambda=-1i*kx*c

        v_kz0=eig_vec_kz0(1:flag.Ny_full-2); %get the v velocity mode.
        u_kz0=-output.D1_bc*v_kz0/(1i*flag.kx_list); %solve u velocity using continuity equation with kz=0

        dpdy_kz0=-(1i*flag.kx_list*(output.U_bar-c_kz0).*v_kz0-1/output.Re*(output.D2_bc-output.kx^2*output.I_bc)*v_kz0);
        p_kz0=(1i*flag.kx_list*(output.U_bar-c_kz0).*u_kz0+v_kz0.*output.d_U_bar-1/output.Re*(output.D2_bc-output.kx^2*output.I_bc)*u_kz0)/(-1i*flag.kx_list);
        
        for y_ind=1:flag.Ny_full-2
            %to check directly performing the integration also match p_kz0
            %computed from x-momentum equation. p_kz0_int is not used. 
            p_kz0_int(y_ind)=sum(dpdy_kz0(y_ind:end).*output.y_diff_unbounded_bc.weight(y_ind:end));
        end

        int_weight=output.y_diff_unbounded_bc.weight; %integration weight.
        Pe=flag.Re_list*flag.Pr_list;
        %below is new for low Prandtl number cases
        for kz_ind=1:length(flag.kz_list)
            kx=flag.kx_list(1);
            kz=flag.kz_list(kz_ind);
            Dp=output.D2_bc*p_kz0-(kx^2+kz^2)*p_kz0;
            %Dp=output.D2_bc*p_kz0-(kx^2)*p_kz0;

            num=kz^2*sum((p_kz0.*Dp)./(output.U_bar-c_kz0).*int_weight);
            den=sum((2*p_kz0.*dpdy_kz0)./(output.U_bar-c_kz0).^4.*output.d_U_bar/(1i*kx).*int_weight);
            %c1_const=flag.kx_list^2/2*sum(p_kz0.^2.*int_weight)/den; %formula of c1 without kz^2
        
            % c1_const_list(Ri_ind)=c1_const; 
            c1_asymptotic(kz_ind)=num/den; %c1 as a function of kz
        
        end
        Fr2=1/flag.Ri_bulk_list(Ri_ind);
        c_asymptotic=c_kz0+Fr2/Pe*c1_asymptotic;
        lambda_asymptotic{Ri_ind}=-1i*flag.kx_list*c_asymptotic;

        % Plot/export raw kz=0 eigenfunctions; retain relative magnitudes.
        f0_kz0=1i*dpdy_kz0/flag.kx_list./((output.U_bar-c_kz0).^2);
        plot_kz0_eigenfunctions(output.y_bc,v_kz0,u_kz0,p_kz0,flag,Ri_ind,flag.kx_list(1), ...
            f0_kz0,dpdy_kz0,output.U_bar,output.d_U_bar,output.dd_U_bar);
    end

    if any(strcmp(flag.post, {'eig_stratified_spanwise_tanh_DCB_figure4a', ...
            'eig_stratified_spanwise_tanh_asymptotic', 'eig_stratified_spanwise_bickley_asymptotic'}))
    
        % Reuse the computed growth rates for each requested horizontal axis.
        plot_axis_names={'kz'};
        if scaled_kz
            plot_axis_names={'kz','Fr_kz'};
        end
        for plot_axis_ind=1:numel(plot_axis_names)
        for plot_component_ind=1:2 % Real and imaginary parts of the same eigenvalues.
                csv_value_name='real_lambda';
                if plot_component_ind==2
                    csv_value_name='imag_lambda';
                end
                csv_axis_name=plot_axis_names{plot_axis_ind};
                plot_kz_list=flag.kz_list;
                data = cell(1,6); % Discard curves from earlier plots.
                data{1}.x=plot_kz_list;
                data{1}.y=output.eig_real_max{1};
                data{2}.x=plot_kz_list;
                data{2}.y=output.eig_real_max{2};
                data{3}.x=plot_kz_list;
                data{3}.y=output.eig_real_max{3};

                %asymptotic prediction based on Eq. (65) 
                data{4}.x=plot_kz_list;
                data{4}.y=real(lambda_asymptotic{1});
                data{5}.x=plot_kz_list;
                data{5}.y=real(lambda_asymptotic{2});
                data{6}.x=plot_kz_list;
                data{6}.y=real(lambda_asymptotic{3});
                if plot_component_ind==2
                    for plot_Ri_ind=1:3
                        data{plot_Ri_ind}.y=output.eig_imag_max{1,1,plot_Ri_ind,1};
                        data{plot_Ri_ind+3}.y=imag(lambda_asymptotic{plot_Ri_ind});
                    end
                end

                if scaled_kz
                    for plot_curve_ind=1:6
                        plot_Fr=1/sqrt(flag.Ri_bulk_list(mod(plot_curve_ind-1,3)+1));
                        if strcmp(csv_axis_name,'Fr_kz')
                            data{plot_curve_ind}.x=flag.Fh_kz_list;
                        else
                            % Each Fr has its own physical kz grid, shared by both methods.
                            data{plot_curve_ind}.x=flag.Fh_kz_list/plot_Fr;
                        end
                    end
                end
                % Discard tanh Fr=1 points beyond Fr*kz=0.82 before plots/CSV export.
                if strcmp(flag.mean,'tanh')
                    for cutoff_curve_ind=1:numel(data)
                        cutoff_Fr=1/sqrt(flag.Ri_bulk_list(mod(cutoff_curve_ind-1,3)+1));
                        if abs(cutoff_Fr-1)<1e-12
                            % For Fr=1, kz and Fr*kz are the same coordinate.
                            cutoff_keep=data{cutoff_curve_ind}.x<=0.82;
                            data{cutoff_curve_ind}.x=data{cutoff_curve_ind}.x(cutoff_keep);
                            data{cutoff_curve_ind}.y=data{cutoff_curve_ind}.y(cutoff_keep);
                        end
                    end
                end
                plot_config.Markerindex=3;
                plot_config.loglog=[0,0];
                plot_config.ylim_list=[1,0,0.2];
                plot_config.ytick_list=[1,0,0.1,0.2];
                plot_config.xlim_list=[1,0,6];
                plot_config.xtick_list=[1,0,1,2,3,4,5,6];
                plot_config.linewidth=2.5;
                plot_config.user_color_style_marker_list = ...
                        {'k-o','r--s','b:^','k-+','r--x','b:d'};
                plot_config.legend_list = {1};
                for plot_curve_ind = 1:6
                    plot_Fr = 1/sqrt(flag.Ri_bulk_list(mod(plot_curve_ind-1,3)+1));
                    plot_curve_types = {'Numerical', 'Asymptotic'};
                    plot_config.legend_list{plot_curve_ind+1} = sprintf('%s, Fr=%g', ...
                        plot_curve_types{ceil(plot_curve_ind/3)}, plot_Fr);
                end
                plot_config.label_list={1,'$k_z$','$\mathrm{Re}[\lambda]$'};
                if scaled_kz
                    if strcmp(csv_axis_name,'Fr_kz')
                        plot_config.label_list={1,'$Fr\,k_z$','$\mathrm{Re}[\lambda]$'};
                        if strcmp(flag.mean,'tanh')
                            plot_config.xlim_list=[1,0,max(flag.Fh_kz_list)];
                        else
                            plot_config.xlim_list=[1,0,4];
                        end
                        plot_config.xtick_list=[1,linspace(0,plot_config.xlim_list(3),5)];
                    else
                        plot_kz_max=max(flag.Fh_kz_list)*sqrt(max(flag.Ri_bulk_list));
                        if strcmp(flag.mean,'tanh')
                            % Leave room to the right of the curves for the legend.
                            plot_kz_max=5*ceil(1.35*plot_kz_max/5);
                        end
                        plot_config.xlim_list=[1,0,plot_kz_max];
                        plot_config.xtick_list=[1,linspace(0,plot_kz_max,5)];
                    end
                end
                plot_config.name=sprintf([flag.mean,'_growth_rate_vs_',csv_axis_name,'_asymptotic_low_Pr_Re_%g_Pe_%g_kx_%g.png'], ...
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), flag.kx_list(1));
                plot_config.name=fullfile(flag.results_folder,plot_config.name);
                % if scaled_kz
                %     % Distinct symbols/styles distinguish even exactly coincident curves.
                %     plot_config.user_color_style_marker_list = ...
                %         {'k-o','r--s','b:^','k-+','r--x','b:d'};
                % end
                if plot_component_ind==2
                    plot_config.label_list{3}='$\mathrm{Im}[\lambda]$';
                    % Imaginary parts may be negative or zero: choose their own limits.
                    plot_imag_values=cellfun(@(curve) curve.y(:),data,'UniformOutput',false);
                    plot_imag_values=vertcat(plot_imag_values{:});
                    plot_imag_values=plot_imag_values(isfinite(plot_imag_values));
                    plot_imag_min=min(plot_imag_values);
                    plot_imag_max=max(plot_imag_values);
                    plot_imag_padding=0.05*max([plot_imag_max-plot_imag_min,abs(plot_imag_min),abs(plot_imag_max),1e-6]);
                    % Reserve additional headroom for the upper-right legend.
                    plot_imag_upper_padding=max(plot_imag_padding,0.4*(plot_imag_max-plot_imag_min));
                    plot_config.ylim_list=[1,plot_imag_min-plot_imag_padding,plot_imag_max+plot_imag_upper_padding];
                    plot_config.ytick_list=[1,linspace(plot_imag_min-plot_imag_padding,plot_imag_max+plot_imag_upper_padding,5)];
                    plot_config.name=strrep(plot_config.name,'.png','_imag.png');
                else
                    plot_config.name=strrep(plot_config.name,'.png','_real.png');
                end
                plot_line(data,plot_config);
                if strcmp(flag.mean,'tanh')
                    comparison_legend=findobj(gcf,'Type','legend');
                    if strcmp(csv_axis_name,'Fr_kz')
                        set(comparison_legend,'BackgroundAlpha',0.5,'Location','southwest');
                    else
                        set(comparison_legend,'BackgroundAlpha',0.5,'Location','northeast');
                    end
                end
                comparison_title = ['$Re=',latex_scientific(flag.Re_list(1)), ...
                    ',\quad Pe=',latex_scientific(flag.Re_list(1)*flag.Pr_list(1)), ...
                    ',\quad k_x=',sprintf('%.4f',flag.kx_list(1)),'$'];
                title(gca,comparison_title,'Interpreter','latex','FontSize',40);
                % if any(strcmp(flag.mean,{'tanh','bickley_jet'}))
                %     % Ten sparse markers per curve, with six interleaved phases.
                %     plot_marker_count = 30;
                %     plot_markers = {'o','s','^','+','x','d'};
                %     plot_handles = flipud(findobj(gca,'Type','line'));
                %     for plot_curve_ind = 1:6
                %         % Match plotted data rather than relying on pre-existing markers.
                %         plot_match = find(arrayfun(@(h) ...
                %             isequal(h.XData(:),data{plot_curve_ind}.x(:)) && ...
                %             isequal(h.YData(:),data{plot_curve_ind}.y(:)),plot_handles),1);
                %         if isempty(plot_match)
                %             continue;
                %         end
                %         plot_handle = plot_handles(plot_match);
                %         plot_handles(plot_match) = []; % Each curve is styled once.
                %         plot_point_count = numel(plot_handle.XData);
                %         plot_marker_phase = (plot_curve_ind-0.5)/6;
                %         plot_marker_fraction = ((0:plot_marker_count-1)+plot_marker_phase)/plot_marker_count;
                %         plot_marker_indices = 1 + round(plot_marker_fraction*(plot_point_count-1));
                %         plot_handle.Marker = plot_markers{plot_curve_ind};
                %         plot_handle.MarkerIndices = unique(plot_marker_indices);
                %         plot_handle.MarkerSize = 10;
                %         plot_handle.MarkerFaceColor = 'none';
                %     end
                %     % Refresh the saved figure after updating markers and legend symbols.
                % end

                % Save the final title and marker settings for either base flow.
                exportgraphics(gcf,plot_config.name,'Resolution',300, ...
                    'Units','pixels','Padding',100,'BackgroundColor','white');

                % Save every plotted curve as paired kz/lambda columns for replotting.
                % Reference curves may have different lengths; unused cells stay blank.
                csv_curve_types = {'Numerical', 'Asymptotic'};
                csv_curve_count = 6;
                csv_row_count = max(cellfun(@(curve) numel(curve.x), data(1:csv_curve_count)));
                csv_export = cell(csv_row_count + 1, 2*csv_curve_count);
                for csv_curve_ind = 1:csv_curve_count
                    csv_Fr = 1/sqrt(flag.Ri_bulk_list(mod(csv_curve_ind-1, 3) + 1));
                    csv_label = sprintf('%s Fr=%g', csv_curve_types{ceil(csv_curve_ind/3)}, csv_Fr);
                    csv_col = 2*csv_curve_ind - 1;
                    csv_export(1, csv_col:csv_col+1) = {[csv_label ', ' csv_axis_name], [csv_label ', ' csv_value_name]};
                    csv_export(2:numel(data{csv_curve_ind}.x)+1, csv_col) = num2cell(data{csv_curve_ind}.x(:));
                    csv_export(2:numel(data{csv_curve_ind}.y)+1, csv_col+1) = num2cell(data{csv_curve_ind}.y(:));
                end
                csv_filename = sprintf([flag.mean,'_growth_rate_vs_',csv_axis_name,'_asymptotic_low_Pr_Re_%g_Pe_%g_kx_%g.csv'], ...
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), flag.kx_list(1));
                if plot_component_ind==2
                    csv_filename=strrep(csv_filename,'.csv','_imag.csv');
                else
                    csv_filename=strrep(csv_filename,'.csv','_real.csv');
                end
                csv_filename=fullfile(flag.results_folder,csv_filename);
                writecell(csv_export, csv_filename);
        end % Real/imaginary comparison versions.
        end % Horizontal-axis versions; no eigenvalues or asymptotics are recomputed.
        
    end 

end
end

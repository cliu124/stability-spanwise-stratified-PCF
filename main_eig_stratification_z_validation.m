clear all;
close all;
clc;

%initialization of different post processing and solution branch
flag.rotation_list=0; %In this code, all rotation number is zero. 

% select one option of flag.post from below. 
% 'eig_stratified_FFLWL_figure5a_validation',...
% 'eig_stratified_FFLWL_figure10_validation',...
% 'eig_stratified_spanwise_LHBLCF_figure6b_validation'
flag.post='eig_stratified_FFLWL_figure5a_validation';

if strcmp(flag.post,'eig_stratified_FFLWL_figure5a_validation')
    
    %%figure 5 of Facchini Favier Le Gal Wang Le Bars 2018 linear_instability_of_the_stratified_plane_couette_flow (1)
    flag.Re_list=[665.86
            700
            977.303
            1172.3
            1424.14
            1739.28
            2789.12
            4010.41
            6077.53
            10019.5
            20074.2
            50102.6
            ];
    flag.Ri_bulk_list=1/0.4^2;
    flag.Pr_list=100000;
    kx_list=linspace(0.01,2,61);
    kz_list=linspace(0.01,30,61);
    Ny_full=82;
    flag.mean='laminar_cou';
elseif strcmp(flag.post,'eig_stratified_FFLWL_figure10_validation')

    % %%-------------------------
    %%figure 10 of Facchini Favier Le Gal Wang Le Bars 2018 linear_instability_of_the_stratified_plane_couette_flow (1)
   flag.Re_list=966;
   flag.Ri_bulk_list=1/0.82^2;
   flag.Pr_list=7;
   kx_list=0.96;
   kz_list=5.16;
   Ny_full=92;
   flag.post_eig_stratified_FFLWL_figure10_validation=1;
    %%---------------------------
   flag.mean='laminar_cou';

elseif strcmp(flag.post,'eig_stratified_spanwise_LHBLCF_figure6b_validation')
    %%figure 6(b) of 
       %%Le Gal P, Harlander U, Borcia ID, Le Dizès S, Chen J, Favier B. Instability of vertically stratified horizontal plane Poiseuille flow. Journal of Fluid Mechanics. 2021 Jan;907.
   flag.Re_list=10000;
   flag.Ri_bulk_list=1/1.1^2;
   flag.Pr_list=1;
   flag.stratified_direction='z';
   kx_list=linspace(0.001,2.6,60);
   kz_list=linspace(0.001,18,60);
   Ny_full=92;
   flag.mean='laminar_poi';

else
    error('Wrong flag.post')
end


[~, DM] = chebdif(Ny_full, 2);
D1_full=DM(:,:,1);
D2_full=DM(:,:,2);

D1_bc=D1_full(2:Ny_full-1,2:Ny_full-1);
D2_bc=D2_full(2:Ny_full-1,2:Ny_full-1);

[y_bc, D4_bc] = cheb4c(Ny_full);

I_bc=eye(Ny_full-2,Ny_full-2);
zero_bc=zeros(Ny_full-2,Ny_full-2);

if strcmp(flag.mean,'laminar_cou')
    U_bar=-y_bc;
    d_U_bar=-ones(Ny_full-2,1);
    dd_U_bar=zeros(Ny_full-2,1);
elseif strcmp(flag.mean,'laminar_poi')
    U_bar=-(1-y_bc.^2);
    d_U_bar=2*y_bc;
    dd_U_bar=2*ones(Ny_full-2,1);
else
    error('Wrong flag.mean');
end

for re_ind=1:length(flag.Re_list)
    Re=flag.Re_list(re_ind);
    for rotation_ind=1:length(flag.rotation_list)
        rotation=flag.rotation_list(rotation_ind);
        for Ri_ind=1:length(flag.Ri_bulk_list)
            Ri_bulk=flag.Ri_bulk_list(Ri_ind);
            for Pr_ind=1:length(flag.Pr_list)
                Pr=flag.Pr_list(Pr_ind);
                for kx_ind=1:length(kx_list)
                    kx=kx_list(kx_ind);
                    for kz_ind=1:length(kz_list)
                        kz=kz_list(kz_ind);
                        K2= kx^2+kz^2; % Total wave number, in convienent for calculation
                        zi=sqrt(-1);
                        
                        %%add the rotating effect here.
                        %%
                        
                        A11=(D4_bc-2*K2*D2_bc+K2^2*I_bc)/Re+diag(dd_U_bar)*zi*kx*I_bc-zi*kx*diag(U_bar)*(D2_bc-K2*I_bc); %%Orr-Sommerfeld operator
                       
                        A21= -zi*kz*diag(d_U_bar)*I_bc; %Coulping operator
                        A22= -zi*kx*diag(U_bar)*I_bc+1/Re*(D2_bc-K2*I_bc); %Squire operator  
                        %%Inverse of [Laplacian, zeros; zeros, I]
                        inv_lap=inv([D2_bc-K2*I_bc, zero_bc; zero_bc, I_bc]);
                        A= inv_lap*[A11, zero_bc-rotation*1i*kz*I_bc; A21+rotation*1i*kz*I_bc, A22];
                        
                        
                        dz_rho_mean=-I_bc; %%For spanwise stratification, I can only work on this version...
                        
                        A_stratified=[A, (-Ri_bulk)*inv_lap*[ -zi*kz*D1_bc; -zi*kx*I_bc];
                                    -dz_rho_mean*[zi*kz*D1_bc, zi*kx*I_bc]/K2,-zi*kx*diag(U_bar)*I_bc+1/Re/Pr*(D2_bc-K2*I_bc)];
                        A=A_stratified;
                        eig_val{re_ind,rotation_ind,Ri_ind,Pr_ind,kx_ind,kz_ind}=eig(A);
                        eig_real_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=max(real(eig(A)));
                            
                    end
                end
                eig_real_max_kx_kz(re_ind,rotation_ind,Ri_ind,Pr_ind)=...
                    max(max(eig_real_max{re_ind,rotation_ind,Ri_ind,Pr_ind}));
            end
        end
    end
end


if strcmp(flag.post,'eig_stratified_FFLWL_figure5a_validation')

%%This is the validation for spanwise stratification
 %%Comparing results of figure 5(a) of 
 %%Facchini G, Favier B, Le Gal P, Wang M, Le Bars M. The linear instability of the stratified plane Couette flow. Journal of Fluid Mechanics. 2018 Oct;853:205-34.
 %%This is for infinite Sc number, plot the growth rate varying over Re
 %%and Froude number
    data{1}.x=flag.Re_list;
    data{1}.y=eig_real_max_kx_kz;
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
    data{1}.x=real(eig_val{1}/(-1i));
    data{1}.y=imag(eig_val{1}/(-1i));
    plot_config.label_list={1,'$Re(\omega)$','$Im(\omega)$'};
    plot_config.xlim_list=[1,-2,2];
    plot_config.ylim_list=[1,-1,0.05];
    plot_config.print_size=[1,1100,900];
    plot_config.loglog=[0,0];
    plot_config.xtick_list=[1,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2];
    plot_config.ytick_list=[1,-1,-0.8,-0.6,-0.4,-0.2,0];
    plot_config.name=['validation_stratification_spanwise_eig_FFLWL_figure10a','.png'];
    plot_config.Markerindex=3;
    plot_config.user_color_style_marker_list={'ksquare','rsquare'};
    plot_line(data,plot_config);
    
    
    data{1}.x=real(eig_val{1}/(-1i));
    data{1}.y=imag(eig_val{1}/(-1i));
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
    plot_config.user_color_style_marker_list={'k*','rsquare'};
    plot_line(data,plot_config);

elseif strcmp(flag.post,'eig_stratified_spanwise_LHBLCF_figure6b_validation')
    data{1}.x=kx_list;
    data{1}.y=kz_list;
    data{1}.z=eig_real_max{1}';
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
    plot_config.user_color_style_marker_list={'k--'};
    plot_contour(data,plot_config);
else
    error('Wrong flag.post')
end

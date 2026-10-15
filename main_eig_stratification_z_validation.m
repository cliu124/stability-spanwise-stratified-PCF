clear all;
close all;
clc;

% Half-transparent backgrounds for every legend, including saved figures.
legend_background_alpha_previous=get(groot,'defaultLegendBackgroundAlpha');
legend_background_cleanup=onCleanup(@() set(groot,'defaultLegendBackgroundAlpha',legend_background_alpha_previous));
set(groot,'defaultLegendBackgroundAlpha',0.5);

%initialization of different post processing and solution branch
flag.rotation_list=0; %In this code, all rotation number is zero. 

% select one option of flag.post from below. 
% 'eig_stratified_FFLWL_figure5a_validation',... Couette flow
% 'eig_stratified_FFLWL_figure10_validation',... Couette flow
% 'eig_stratified_spanwise_LHBLCF_figure6b_validation',... Poiseuille flow
% 'eig_stratified_spanwise_tanh_DCB_figure4a',... tanh profile
% 'eig_stratified_spanwise_bickley_DCB_figure5'% Bickley Jet

flag.asymptotic=0;
flag.asymptotic_low_Pr=0;

%flag.post='eig_stratified_spanwise_tanh_DCB_figure4a'; 
%flag.post='eig_stratified_spanwise_bickley_DCB_figure5';
%flag.post='eig_stratified_spanwise_tanh_growth_rate_kx';

flag.post='eig_stratified_spanwise_tanh_asymptotic';
%flag.post='eig_stratified_spanwise_bickley_asymptotic'; 
%flag.post='eig_stratified_spanwise_bickley_growth_rate_kx';
%flag.post='eig_stratified_spanwise_bickley_CLDM_figure1';

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
   
elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_DCB_figure4a')
   flag.Re_list=[100000];
   flag.Pr_list=[0.000001];
   flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
   kx_list=0.4449;
   if flag.Pr_list<0.01
       kz_list=linspace(0,5,100);
   else
       kz_list=linspace(0,18,100);
   end
   Ny_full=202; %%32.
   flag.mean='tanh';

elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_DCB_figure5')
    flag.Re_list=[100000];
    flag.Pr_list=[1];
    flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
    kx_list=0.9021;%linspace(0.01,2,20);
    Fh_kz_list=linspace(0.01,4,20);
    % 
    % flag.Ri_bulk_list=1;
    % kx_list=0.9023;
    % Fh_kz_list=0;

    Ny_full=402; 
    flag.mean='bickley_jet';
elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_CLDM_figure1')
    flag.Re_list=100000;
    flag.Pr_list=1;
    flag.Ri_bulk_list=1;
    kx_list=linspace(0.01,2,60);
    Fh_kz_list=linspace(0.01,4,60);
    Ny_full=202;
    flag.mean='bickley_jet';
elseif strcmp(flag.post, 'eig_stratified_spanwise_tanh_growth_rate_kx')
    flag.Re_list=[100000];
    flag.Pr_list=[1];
    flag.Ri_bulk_list=1./[0.1,0.01].^2;
    kx_list=linspace(0.001,1.05,200);
    kz_list=[0,2,4,6,8];
    Ny_full=202; %%32.
    flag.mean='tanh';
elseif strcmp(flag.post, 'eig_stratified_spanwise_bickley_growth_rate_kx')
    flag.Re_list=1e8;
    flag.Pr_list=[1];
    flag.Ri_bulk_list=1./[0.1,0.01].^2;
    kx_list=linspace(0.001,2,200);
    %kx_list=linspace(1.9,2,10);
    
    kz_list=[0,2,4,6,8];
    Ny_full=302; %%32.
    flag.mean='bickley_jet';

elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_asymptotic')
    %compare the asymptotic with the numerical results for tanh profile
    flag.Re_list=[1e5];
    %flag.Pr_list=[0.000001];
    flag.Pr_list=1;
    %flag.Pr_list=1e-5;
    %flag.Pr_list=1e-6;
    flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
    kx_list=0.4449;
    if flag.Pr_list<0.01
        kz_list=linspace(0,5,100);
        flag.asymptotic_low_Pr=1;
        flag.asymptotic=0;
    else
        kz_list=linspace(0,18,100);
        %kz_list=18;
        %kz_list=kz_list(10);
        flag.asymptotic_low_Pr=0;
        flag.asymptotic=1;
    end
    Ny_full=203; %%32.
    flag.mean='tanh';
elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
    %compare the asymptotic with the numerical results for bickley jet
    %profile. 
    flag.Re_list=[1e8];
    %flag.Pr_list=[0.000001];
    flag.Pr_list=1;
    %flag.Pr_list=1e-8;
    %flag.Pr_list=1e-9;
    flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
    kx_list=0.9021;
    if flag.Pr_list<0.01
        flag.asymptotic_low_Pr=1;
        flag.asymptotic=0;
        Fh_kz_list=linspace(0,2,50); % Common horizontal coordinate Fr*kz.

    else
        Fh_kz_list=linspace(0,4,100); % Common horizontal coordinate Fr*kz.
        flag.asymptotic_low_Pr=0;
        flag.asymptotic=1;
    end
    kz_list=Fh_kz_list; % Converted to physical kz separately for each Fr.

    Ny_full=302; %%32.
    flag.mean='bickley_jet';

else
    error('Wrong flag.post')
end


% Group outputs for each single mean-flow/Re/Pe parameter combination.
if isscalar(flag.Re_list) && isscalar(flag.Pr_list)
    flag.results_folder=fullfile(pwd,sprintf('%s_Re_%g_Pe_%g', ...
        flag.mean,flag.Re_list,flag.Re_list*flag.Pr_list));
    if ~isfolder(flag.results_folder)
        [created,message]=mkdir(flag.results_folder);
        if ~created
            error('Could not create results folder %s: %s',flag.results_folder,message);
        end
    end
else
    % Existing validation plots over parameter sweeps combine several cases.
    flag.results_folder=pwd;
end

[eig_val,eig_vec,eig_real_max,eig_imag_max,eig_real_max_kx_kz]=eig_solve(flag,Ny_full,kx_list,kz_list)


function [eig_val,eig_vec,eig_real_max,eig_imag_max,eig_real_max_kx_kz]=eig_solve(flag,Ny_full,kx_list,kz_list)

[y_list_full, DM] = chebdif(Ny_full, 4);
D1_full=DM(:,:,1);
D2_full=DM(:,:,2);
D3_full=DM(:,:,3);
D4_full=DM(:,:,4);

%construct y_diff_full that will be used for tanh and bickley jet
y_diff_full.D1=D1_full;
y_diff_full.D2=D2_full;
y_diff_full.D3=D3_full;
y_diff_full.D4=D4_full;
y_diff_full.y_list=y_list_full;


[~,w]=clencurt(Ny_full-1);
y_diff_full.weight=w';
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
elseif strcmp(flag.mean,'tanh') 
    %conduct the coordinate transform suggested by 
    %Deloncle A, Chomaz JM, Billant P. Three-dimensional stability of a horizontally sheared flow in a stably stratified fluid. Journal of Fluid Mechanics. 2007 Jan;570:297-305.
    [y_diff_unbounded_full,y_diff_unbounded_bc]=channel2unbounded(y_diff_full);
    D1_bc=y_diff_unbounded_bc.D1;
    D2_bc=y_diff_unbounded_bc.D2;
    D4_bc=y_diff_unbounded_bc.D4;
    
    y_bc=y_diff_unbounded_bc.y_list;
    syms x 
    U_sym=tanh(x);
    dU_sym=diff(U_sym);
    ddU_sym=diff(dU_sym);

    U_bar=double(subs(U_sym,y_bc));
    d_U_bar=double(subs(dU_sym,y_bc));
    dd_U_bar=double(subs(ddU_sym,y_bc));      

elseif strcmp(flag.mean,'bickley_jet')    
    %conduct the coordinate transform suggested by 
    %Deloncle A, Chomaz JM, Billant P. Three-dimensional stability of a horizontally sheared flow in a stably stratified fluid. Journal of Fluid Mechanics. 2007 Jan;570:297-305.
    [y_diff_unbounded_full,y_diff_unbounded_bc]=channel2unbounded(y_diff_full);
    D1_bc=y_diff_unbounded_bc.D1;
    D2_bc=y_diff_unbounded_bc.D2;
    D4_bc=y_diff_unbounded_bc.D4;
    y_bc=y_diff_unbounded_bc.y_list;
    syms x 
    U_sym=sech(x)^2;
    dU_sym=diff(U_sym);
    ddU_sym=diff(dU_sym);

    U_bar=double(subs(U_sym,y_bc));
    d_U_bar=double(subs(dU_sym,y_bc));
    dd_U_bar=double(subs(ddU_sym,y_bc));         
else
    error('Wrong flag.mean');
end

for re_ind=1:length(flag.Re_list)
    Re=flag.Re_list(re_ind);
    for rotation_ind=1:length(flag.rotation_list)
        rotation=flag.rotation_list(rotation_ind);
        for Ri_ind=1:length(flag.Ri_bulk_list)
            Ri_bulk=flag.Ri_bulk_list(Ri_ind);
            if strcmp(flag.mean,'bickley_jet') && ~strcmp(flag.post,'eig_stratified_spanwise_bickley_growth_rate_kx')
                Fh=1/sqrt(Ri_bulk);
                kz_list=Fh_kz_list/Fh;
            end
            
            for Pr_ind=1:length(flag.Pr_list)
                Pr=flag.Pr_list(Pr_ind);
                for kx_ind=1:length(kx_list)
                    kx=kx_list(kx_ind);
                    for kz_ind=1:length(kz_list)
                        kz=kz_list(kz_ind);
                        K2= kx^2+kz^2; % Total wave number, in convienent for calculation
                        zi=sqrt(-1);
                        
                        
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
                        [V,D]=eig(A);
                        eig_val{re_ind,rotation_ind,Ri_ind,Pr_ind,kx_ind,kz_ind}=diag(D);
                     
                        [eig_real_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind),ind]=max(real(diag(D)));                          
                        %eig_imag_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=imag(D(ind,ind));
                        % 
                        if max(real(diag(D)))>0
                            eig_imag_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=imag(D(ind,ind));
                        else 
                            eig_imag_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=NaN;
                        end
                        eig_vec{re_ind,rotation_ind,Ri_ind,Pr_ind,kx_ind,kz_ind}=V(:,ind);
                    end
                end
                eig_real_max_kx_kz(re_ind,rotation_ind,Ri_ind,Pr_ind)=...
                    max(max(eig_real_max{re_ind,rotation_ind,Ri_ind,Pr_ind}));
            end
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
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
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
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
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
    plot_config.name=fullfile(flag.results_folder,plot_config.name);
    plot_config.user_color_style_marker_list={'k--'};
    plot_contour(data,plot_config);
    
    
elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_DCB_figure4a') || ...
        (strcmp(flag.post,'eig_stratified_spanwise_tanh_asymptotic') && strcmp(flag.mean,'tanh'))
    DCB_figure4a_Fh1=[0	0.189744
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
    DCB_figure4a_Fh01=[0	0.189714
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
    
    DCB_figure4a_Fh005=[0	0.189714
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
    
    data{1}.x=kz_list;
    data{1}.y=eig_real_max{1};
    data{2}.x=kz_list;
    data{2}.y=eig_real_max{2};
    data{3}.x=kz_list;
    data{3}.y=eig_real_max{3};
    data{4}.x=DCB_figure4a_Fh1(:,1);
    data{4}.y=DCB_figure4a_Fh1(:,2);
    data{5}.x=DCB_figure4a_Fh01(:,1);
    data{5}.y=DCB_figure4a_Fh01(:,2);
    data{6}.x=DCB_figure4a_Fh005(:,1);
    data{6}.y=DCB_figure4a_Fh005(:,2);
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
    data{1}.x=Fh_kz_list;%kz_list/sqrt(flag.Ri_bulk_list(1));
    data{1}.y=eig_real_max{1};
    data{2}.x=Fh_kz_list;%kz_list/sqrt(flag.Ri_bulk_list(2));
    data{2}.y=eig_real_max{2};
    data{3}.x=Fh_kz_list;%kz_list/sqrt(flag.Ri_bulk_list(3));
    data{3}.y=eig_real_max{3};
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
        for kz_ind=1:length(kz_list)
            eig_real_max{Ri_ind};
            data{kz_ind}.x=kx_list;
            data{kz_ind}.y=eig_real_max{Ri_ind}(:,kz_ind);
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
        csv_curve_count = length(kz_list);
        csv_export = cell(numel(kx_list) + 1, 2*csv_curve_count);
        for csv_curve_ind = 1:csv_curve_count
            csv_label = sprintf('Numerical Fr=%g, kz=%g', Fr, kz_list(csv_curve_ind));
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
    growth_rate = real(eig_real_max{Ri_ind}).';

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

    [C, h] = contour(kx_list, kz_list, growth_rate, levels, ...
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

% save('data_tanh');

if flag.asymptotic
    for Ri_ind=1:length(flag.Ri_bulk_list)
        if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
            % Evaluate asymptotics on the same physical kz grid as the eigensolver.
            kz_list=Fh_kz_list*sqrt(flag.Ri_bulk_list(Ri_ind));
        end

        %Get the eigenvector associated with kz=0. Here assume the first
        %value in kz_list is 0. Thus just take index=1 for the kz. 
        eig_vec_kz0=eig_vec{1,1,Ri_ind,1,1,1};
        eig_val_kz0=eig_val{1,1,Ri_ind,1,1,1};
        [val,ind]=max(real(eig_val_kz0)); %get the index with the maximal real part
        
        lambda_kz0=eig_val_kz0(ind); %get the lambda for most unstable mode
        c_kz0=lambda_kz0/(-1i*kx_list); %get corresponding c valie, lambda=-1i*kx*c

        v_kz0=eig_vec_kz0(1:Ny_full-2); %get the v velocity mode. 
        u_kz0=-D1_bc*v_kz0/(1i*kx_list); %solve u velocity using continuity equation with kz=0

        dpdy_kz0=-(1i*kx_list*(U_bar-c_kz0).*v_kz0-1/Re*(D2_bc-kx^2*I_bc)*v_kz0);
        p_kz0=(1i*kx_list*(U_bar-c_kz0).*u_kz0+v_kz0.*d_U_bar-1/Re*(D2_bc-kx^2*I_bc)*u_kz0)/(-1i*kx_list);
        
        f0_kz0=1i*dpdy_kz0/kx_list./((U_bar-c_kz0).^2);
        
        for y_ind=1:Ny_full-2
            %to check directly performing the integration also match p_kz0
            %computed from x-momentum equation. p_kz0_int is not used. 
            p_kz0_int(y_ind)=sum(dpdy_kz0(y_ind:end).*y_diff_unbounded_bc.weight(y_ind:end));
        end

        int_weight=y_diff_unbounded_bc.weight; %integration weight. 
        den=sum((dpdy_kz0.^2+kx_list^2*p_kz0.^2)./(U_bar-c_kz0).^3.*int_weight);
        c1_const=kx_list^2/2*sum(p_kz0.^2.*int_weight)/den; %formula of c1 without kz^2
        
        c1_const_list(Ri_ind)=c1_const; 

        c1_asymptotic=c1_const*kz_list.^2; %c1 as a function of kz
        Fr2=1/flag.Ri_bulk_list(Ri_ind);
        c_asymptotic=c_kz0+Fr2*c1_asymptotic;
        lambda_asymptotic{Ri_ind}=-1i*kx_list*c_asymptotic;

        % Plot/export raw kz=0 eigenfunctions; retain relative magnitudes.
        plot_kz0_eigenfunctions(y_bc,v_kz0,u_kz0,p_kz0,flag,Ri_ind,kx_list(1), ...
            f0_kz0,dpdy_kz0,U_bar,d_U_bar,dd_U_bar);
    end

    if any(strcmp(flag.post, {'eig_stratified_spanwise_tanh_DCB_figure4a', ...
            'eig_stratified_spanwise_tanh_asymptotic', 'eig_stratified_spanwise_bickley_asymptotic'}))
    
        % Reuse the computed growth rates for each requested horizontal axis.
        plot_axis_names={'kz'};
        if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
            plot_axis_names={'kz','Fr_kz'};
        end
        for plot_axis_ind=1:numel(plot_axis_names)
        for plot_component_ind=1:2 % Real and imaginary parts of the same eigenvalues.
                csv_value_name='real_lambda';
                if plot_component_ind==2
                    csv_value_name='imag_lambda';
                end
                csv_axis_name=plot_axis_names{plot_axis_ind};
                plot_kz_list=kz_list;
                data = cell(1,6); % Discard curves from earlier plots.
                data{1}.x=plot_kz_list;
                data{1}.y=eig_real_max{1};
                data{2}.x=plot_kz_list;
                data{2}.y=eig_real_max{2};
                data{3}.x=plot_kz_list;
                data{3}.y=eig_real_max{3};
                % Asymptotic curves precede reference markers in plots and CSV files.
                data{4}.x=plot_kz_list;
                data{4}.y=real(lambda_asymptotic{1});
                data{5}.x=plot_kz_list;
                data{5}.y=real(lambda_asymptotic{2});
                data{6}.x=plot_kz_list;
                data{6}.y=real(lambda_asymptotic{3});
                if plot_component_ind==2
                    for plot_Ri_ind=1:3
                        data{plot_Ri_ind}.y=eig_imag_max{1,1,plot_Ri_ind,1};
                        data{plot_Ri_ind+3}.y=imag(lambda_asymptotic{plot_Ri_ind});
                    end
                end

                if strcmp(flag.mean,'tanh') && plot_component_ind==1
                    data{7}.x=DCB_figure4a_Fh1(:,1);
                    data{7}.y=DCB_figure4a_Fh1(:,2);
                    data{8}.x=DCB_figure4a_Fh01(:,1);
                    data{8}.y=DCB_figure4a_Fh01(:,2);
                    data{9}.x=DCB_figure4a_Fh005(:,1);
                    data{9}.y=DCB_figure4a_Fh005(:,2);
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

                if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
                    for plot_curve_ind=1:6
                        plot_Fr=1/sqrt(flag.Ri_bulk_list(mod(plot_curve_ind-1,3)+1));
                        if strcmp(csv_axis_name,'Fr_kz')
                            data{plot_curve_ind}.x=Fh_kz_list;
                        else
                            % Each Fr has its own physical kz grid, shared by both methods.
                            data{plot_curve_ind}.x=Fh_kz_list/plot_Fr;
                        end
                    end
                end
                plot_config.Markerindex=3;
                plot_config.loglog=[0,0];
                plot_config.ylim_list=[1,0.02,0.2];
                plot_config.ytick_list=[1,0,0.1,0.2];
                plot_config.xlim_list=[1,0,25];
                plot_config.xtick_list=[1,0,5,10,15,20,25];
                plot_config.linewidth=2.5;
                plot_config.fontsize_legend=32;
                plot_config.label_list={1,'$k_z$','$\mathrm{Re}[\lambda]$'};
                if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
                    if strcmp(csv_axis_name,'Fr_kz')
                        plot_config.label_list={1,'$Fr\,k_z$','$\mathrm{Re}[\lambda]$'};
                        plot_config.xlim_list=[1,0,4];
                        plot_config.xtick_list=[1,0,1,2,3,4];
                    else
                        plot_kz_max=max(Fh_kz_list)*sqrt(max(flag.Ri_bulk_list));
                        plot_config.xlim_list=[1,0,plot_kz_max];
                        plot_config.xtick_list=[1,linspace(0,plot_kz_max,5)];
                    end
                end
                plot_config.name= sprintf([flag.mean,'_growth_rate_vs_',csv_axis_name,'_asymptotic_Re_%g_Pe_%g_kx_%g.png'], ...
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), kx_list(1));
                plot_config.name=fullfile(flag.results_folder,plot_config.name);
                % if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
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
                comparison_title = ['$Re=',latex_scientific(flag.Re_list(1)), ...
                    ',\quad Pe=',latex_scientific(flag.Re_list(1)*flag.Pr_list(1)), ...
                    ',\quad k_x=',sprintf('%.4f',kx_list(1)),'$'];
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
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), kx_list(1));
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
        if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
            % Evaluate asymptotics on the same physical kz grid as the eigensolver.
            kz_list=Fh_kz_list*sqrt(flag.Ri_bulk_list(Ri_ind));
        end

        %Get the eigenvector associated with kz=0. Here assume the first
        %value in kz_list is 0. Thus just take index=1 for the kz. 
        eig_vec_kz0=eig_vec{1,1,Ri_ind,1,1,1};
        eig_val_kz0=eig_val{1,1,Ri_ind,1,1,1};
        [val,ind]=max(real(eig_val_kz0)); %get the index with the maximal real part
        
        lambda_kz0=eig_val_kz0(ind); %get the lambda for most unstable mode
        c_kz0=lambda_kz0/(-1i*kx_list); %get corresponding c valie, lambda=-1i*kx*c

        v_kz0=eig_vec_kz0(1:Ny_full-2); %get the v velocity mode. 
        u_kz0=-D1_bc*v_kz0/(1i*kx_list); %solve u velocity using continuity equation with kz=0

        dpdy_kz0=-(1i*kx_list*(U_bar-c_kz0).*v_kz0-1/Re*(D2_bc-kx^2*I_bc)*v_kz0);
        p_kz0=(1i*kx_list*(U_bar-c_kz0).*u_kz0+v_kz0.*d_U_bar-1/Re*(D2_bc-kx^2*I_bc)*u_kz0)/(-1i*kx_list);
        
        for y_ind=1:Ny_full-2
            %to check directly performing the integration also match p_kz0
            %computed from x-momentum equation. p_kz0_int is not used. 
            p_kz0_int(y_ind)=sum(dpdy_kz0(y_ind:end).*y_diff_unbounded_bc.weight(y_ind:end));
        end

        int_weight=y_diff_unbounded_bc.weight; %integration weight. 
        Pe=flag.Re_list*flag.Pr_list;
        %below is new for low Prandtl number cases
        for kz_ind=1:length(kz_list)
            kx=kx_list(1);
            kz=kz_list(kz_ind);
            Dp=D2_bc*p_kz0-(kx^2+kz^2)*p_kz0;
            %Dp=D2_bc*p_kz0-(kx^2)*p_kz0;

            num=kz^2*sum((p_kz0.*Dp)./(U_bar-c_kz0).*int_weight);
            den=sum((2*p_kz0.*dpdy_kz0)./(U_bar-c_kz0).^4.*d_U_bar/(1i*kx).*int_weight);
            %c1_const=kx_list^2/2*sum(p_kz0.^2.*int_weight)/den; %formula of c1 without kz^2
        
            % c1_const_list(Ri_ind)=c1_const; 
            c1_asymptotic(kz_ind)=num/den; %c1 as a function of kz
        
        end
        Fr2=1/flag.Ri_bulk_list(Ri_ind);
        c_asymptotic=c_kz0+Fr2/Pe*c1_asymptotic;
        lambda_asymptotic{Ri_ind}=-1i*kx_list*c_asymptotic;

        % Plot/export raw kz=0 eigenfunctions; retain relative magnitudes.
        f0_kz0=1i*dpdy_kz0/kx_list./((U_bar-c_kz0).^2);
        plot_kz0_eigenfunctions(y_bc,v_kz0,u_kz0,p_kz0,flag,Ri_ind,kx_list(1), ...
            f0_kz0,dpdy_kz0,U_bar,d_U_bar,dd_U_bar);
    end

    if any(strcmp(flag.post, {'eig_stratified_spanwise_tanh_DCB_figure4a', ...
            'eig_stratified_spanwise_tanh_asymptotic', 'eig_stratified_spanwise_bickley_asymptotic'}))
    
        % Reuse the computed growth rates for each requested horizontal axis.
        plot_axis_names={'kz'};
        if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
            plot_axis_names={'kz','Fr_kz'};
        end
        for plot_axis_ind=1:numel(plot_axis_names)
        for plot_component_ind=1:2 % Real and imaginary parts of the same eigenvalues.
                csv_value_name='real_lambda';
                if plot_component_ind==2
                    csv_value_name='imag_lambda';
                end
                csv_axis_name=plot_axis_names{plot_axis_ind};
                plot_kz_list=kz_list;
                data = cell(1,6); % Discard curves from earlier plots.
                data{1}.x=plot_kz_list;
                data{1}.y=eig_real_max{1};
                data{2}.x=plot_kz_list;
                data{2}.y=eig_real_max{2};
                data{3}.x=plot_kz_list;
                data{3}.y=eig_real_max{3};

                %asymptotic prediction based on Eq. (65) 
                data{4}.x=plot_kz_list;
                data{4}.y=real(lambda_asymptotic{1});
                data{5}.x=plot_kz_list;
                data{5}.y=real(lambda_asymptotic{2});
                data{6}.x=plot_kz_list;
                data{6}.y=real(lambda_asymptotic{3});
                if plot_component_ind==2
                    for plot_Ri_ind=1:3
                        data{plot_Ri_ind}.y=eig_imag_max{1,1,plot_Ri_ind,1};
                        data{plot_Ri_ind+3}.y=imag(lambda_asymptotic{plot_Ri_ind});
                    end
                end

                if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
                    for plot_curve_ind=1:6
                        plot_Fr=1/sqrt(flag.Ri_bulk_list(mod(plot_curve_ind-1,3)+1));
                        if strcmp(csv_axis_name,'Fr_kz')
                            data{plot_curve_ind}.x=Fh_kz_list;
                        else
                            % Each Fr has its own physical kz grid, shared by both methods.
                            data{plot_curve_ind}.x=Fh_kz_list/plot_Fr;
                        end
                    end
                end
                plot_config.Markerindex=3;
                plot_config.loglog=[0,0];
                plot_config.ylim_list=[1,0.02,0.2];
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
                if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
                    if strcmp(csv_axis_name,'Fr_kz')
                        plot_config.label_list={1,'$Fr\,k_z$','$\mathrm{Re}[\lambda]$'};
                        plot_config.xlim_list=[1,0,4];
                        plot_config.xtick_list=[1,0,1,2,3,4];
                    else
                        plot_kz_max=max(Fh_kz_list)*sqrt(max(flag.Ri_bulk_list));
                        plot_config.xlim_list=[1,0,plot_kz_max];
                        plot_config.xtick_list=[1,linspace(0,plot_kz_max,5)];
                    end
                end
                plot_config.name=sprintf([flag.mean,'_growth_rate_vs_',csv_axis_name,'_asymptotic_low_Pr_Re_%g_Pe_%g_kx_%g.png'], ...
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), kx_list(1));
                plot_config.name=fullfile(flag.results_folder,plot_config.name);
                % if strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic')
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
                comparison_title = ['$Re=',latex_scientific(flag.Re_list(1)), ...
                    ',\quad Pe=',latex_scientific(flag.Re_list(1)*flag.Pr_list(1)), ...
                    ',\quad k_x=',sprintf('%.4f',kx_list(1)),'$'];
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
                    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1), kx_list(1));
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

% Restore the graphics default after this script finishes.
clear legend_background_cleanup;

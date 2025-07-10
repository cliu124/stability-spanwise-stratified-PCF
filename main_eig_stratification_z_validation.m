clear all;
close all;
clc;

%initialization of different post processing and solution branch
flag.rotation_list=0; %In this code, all rotation number is zero. 

% select one option of flag.post from below. 
% 'eig_stratified_FFLWL_figure5a_validation',... Couette flow
% 'eig_stratified_FFLWL_figure10_validation',... Couette flow
% 'eig_stratified_spanwise_LHBLCF_figure6b_validation',... Poiseuille flow
% 'eig_stratified_spanwise_tanh_DCB_figure4a',... tanh profile
% 'eig_stratified_spanwise_bickley_DCB_figure5'% Bickley Jet

flag.post='eig_stratified_spanwise_tanh_DCB_figure4a'; 
% flag.post='eig_stratified_spanwise_bickley_DCB_figure5';

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
   flag.Pr_list=[1];
   flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
   kx_list=0.4449;
   kz_list=linspace(0,18,100);
   Ny_full=202; %%32.
   flag.mean='tanh';
elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_DCB_figure5')
    flag.Re_list=[10000000];
    flag.Pr_list=[1];
    flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
    kx_list=0.9021;%linspace(0.01,2,20);
    Fh_kz_list=linspace(0.01,4,20);
    Ny_full=402; 
    flag.mean='bickley_jet';
   
else
    error('Wrong flag.post')
end


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
            if strcmp(flag.mean,'bickley_jet') 
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
    
    
elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_DCB_figure4a')
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
    plot_config.label_list={1,'$F_h k_z$','$\sigma$'};
    plot_config.name=['eig_stratified_spanwise_validation_stratification_z_tanh_bickley_DCB_figure5','.png'];
    plot_line(data,plot_config);

else
    error('Wrong flag.post')
end

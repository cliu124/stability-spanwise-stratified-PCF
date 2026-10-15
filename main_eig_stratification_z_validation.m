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
flag.post='eig_stratified_spanwise_tanh_growth_rate_kx';

%flag.post='eig_stratified_spanwise_tanh_asymptotic';
flag.solve=1;
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
    flag.kx_list=linspace(0.01,2,61);
    flag.kz_list=linspace(0.01,30,61);
    flag.Ny_full=82;
    flag.mean='laminar_cou';
elseif strcmp(flag.post,'eig_stratified_FFLWL_figure10_validation')

    % %%-------------------------
    %%figure 10 of Facchini Favier Le Gal Wang Le Bars 2018 linear_instability_of_the_stratified_plane_couette_flow (1)
   flag.Re_list=966;
   flag.Ri_bulk_list=1/0.82^2;
   flag.Pr_list=7;
   flag.kx_list=0.96;
   flag.kz_list=5.16;
   flag.Ny_full=92;
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
   flag.kx_list=linspace(0.001,2.6,60);
   flag.kz_list=linspace(0.001,18,60);
   flag.Ny_full=92;
   flag.mean='laminar_poi';
   
elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_DCB_figure4a')
   flag.Re_list=[100000];
   flag.Pr_list=[0.000001];
   flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
   flag.kx_list=0.4449;
   if flag.Pr_list<0.01
       flag.kz_list=linspace(0,5,100);
   else
       flag.kz_list=linspace(0,18,100);
   end
   flag.Ny_full=202; %%32.
   flag.mean='tanh';

elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_DCB_figure5')
    flag.Re_list=[100000];
    flag.Pr_list=[1];
    flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
    flag.kx_list=0.9021;%linspace(0.01,2,20);
    flag.Fh_kz_list=linspace(0.01,4,20);
    %
    % flag.Ri_bulk_list=1;
    % flag.kx_list=0.9023;
    % flag.Fh_kz_list=0;

    flag.Ny_full=402;
    flag.mean='bickley_jet';
elseif strcmp(flag.post,'eig_stratified_spanwise_bickley_CLDM_figure1')
    flag.Re_list=100000;
    flag.Pr_list=1;
    flag.Ri_bulk_list=1;
    flag.kx_list=linspace(0.01,2,60);
    flag.Fh_kz_list=linspace(0.01,4,60);
    flag.Ny_full=202;
    flag.mean='bickley_jet';
elseif strcmp(flag.post, 'eig_stratified_spanwise_tanh_growth_rate_kx')
    flag.Re_list=[100000];
    flag.Pr_list=[1];
    flag.Ri_bulk_list=1./[0.1,0.01].^2;
    flag.kx_list=linspace(0.001,1.05,200);
    flag.kz_list=[0,2,4,6,8];
    flag.Ny_full=302; %%32.
    flag.mean='tanh';
elseif strcmp(flag.post, 'eig_stratified_spanwise_bickley_growth_rate_kx')
    flag.Re_list=1e8;
    flag.Pr_list=[1];
    flag.Ri_bulk_list=1./[0.1,0.01].^2;
    flag.kx_list=linspace(0.001,2,200);
    %flag.kx_list=linspace(1.9,2,10);
    
    flag.kz_list=[0,2,4,6,8];
    flag.Ny_full=302; %%32.
    flag.mean='bickley_jet';

elseif strcmp(flag.post,'eig_stratified_spanwise_tanh_asymptotic')
    %compare the asymptotic with the numerical results for tanh profile
    flag.Re_list=[1e5];
    %flag.Pr_list=[0.000001];
    flag.Pr_list=1;
    %flag.Pr_list=1e-5;
    %flag.Pr_list=1e-6;
    flag.Ri_bulk_list=1./[1,0.1,0.05].^2;
    flag.kx_list=0.4449;
    if flag.Pr_list<0.01
        %flag.kz_list=linspace(0,5,100);
        flag.Fh_kz_list=linspace(0,0.9,100); % Common horizontal coordinate Fr*kz.
        flag.asymptotic_low_Pr=1;
        flag.asymptotic=0;
    else
        % flag.kz_list=linspace(0,18,100);
        flag.Fh_kz_list=linspace(0,0.9,100); % Common horizontal coordinate Fr*kz.
        %flag.kz_list=18;
        %flag.kz_list=flag.kz_list(10);
        flag.asymptotic_low_Pr=0;
        flag.asymptotic=1;
    end
    flag.kz_list=flag.Fh_kz_list; % Converted to physical kz separately for each Fr.

    flag.Ny_full=302; %%32.
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
    flag.kx_list=0.9021;
    if flag.Pr_list<0.01
        flag.asymptotic_low_Pr=1;
        flag.asymptotic=0;
        flag.Fh_kz_list=linspace(0,2,50); % Common horizontal coordinate Fr*kz.
    else
        flag.Fh_kz_list=linspace(0,4,100); % Common horizontal coordinate Fr*kz.
        flag.asymptotic_low_Pr=0;
        flag.asymptotic=1;
    end
    flag.kz_list=flag.Fh_kz_list; % Converted to physical kz separately for each Fr.

    flag.Ny_full=302; %%32.
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

dataset_file=sprintf([flag.post,'_eig_Re_%g_Pe_%g.mat'], ...
    flag.Re_list(1), flag.Re_list(1)*flag.Pr_list(1));

if flag.solve
    flag.store_all_eig_vec=false;
    output=eig_solve(flag);
    % flag_high_Ny=flag;
    % flag_high_Ny.Ny_full=flag.Ny_full+70;
    % flag_high_Ny.store_all_eig_vec=false;
    % output_high_Ny=eig_solve(flag_high_Ny);
    %
    % output_no_spurious=filter_spurious_modes(output,output_high_Ny,0.01);
    save(dataset_file,'output','flag','-v7.3');

    %save(dataset_file,'output','output_high_Ny','output_no_spurious','flag','flag_high_Ny','-v7.3');
else
    saved_dataset=load(dataset_file,'output','output_high_Ny','output_no_spurious');
    output=saved_dataset.output;
    % if isfield(saved_dataset,'output_high_Ny')
    %     output_high_Ny=saved_dataset.output_high_Ny;
    % end
    % if isfield(saved_dataset,'output_no_spurious')
    %     output_no_spurious=saved_dataset.output_no_spurious;
    % end
    %flag=saved_dataset.flag;
end
% Preserve the final physical grid used by the solver for postprocessing.
flag.kz_list=output.kz_list;

post_validation(flag,output);

post_asymptotics(flag,output);

% Restore the graphics default after this script finishes.
clear legend_background_cleanup;

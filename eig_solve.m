function output=eig_solve(flag)
% Solve the configured stability problem and return spectra and base-flow data.

if ~isfield(flag,'kz_list')
    flag.kz_list=flag.Fh_kz_list;
end

[y_list_full, DM] = chebdif(flag.Ny_full, 4);
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


[~,w]=clencurt(flag.Ny_full-1);
y_diff_full.weight=w';
output.D1_bc=D1_full(2:flag.Ny_full-1,2:flag.Ny_full-1);
output.D2_bc=D2_full(2:flag.Ny_full-1,2:flag.Ny_full-1);

[output.y_bc, output.D4_bc] = cheb4c(flag.Ny_full);

output.I_bc=eye(flag.Ny_full-2,flag.Ny_full-2);
zero_bc=zeros(flag.Ny_full-2,flag.Ny_full-2);

if strcmp(flag.mean,'laminar_cou')
    output.U_bar=-output.y_bc;
    output.d_U_bar=-ones(flag.Ny_full-2,1);
    output.dd_U_bar=zeros(flag.Ny_full-2,1);
elseif strcmp(flag.mean,'laminar_poi')
    output.U_bar=-(1-output.y_bc.^2);
    output.d_U_bar=2*output.y_bc;
    output.dd_U_bar=2*ones(flag.Ny_full-2,1);
elseif strcmp(flag.mean,'tanh') 
    %conduct the coordinate transform suggested by 
    %Deloncle A, Chomaz JM, Billant P. Three-dimensional stability of a horizontally sheared flow in a stably stratified fluid. Journal of Fluid Mechanics. 2007 Jan;570:297-305.
    [y_diff_unbounded_full,output.y_diff_unbounded_bc]=channel2unbounded(y_diff_full);
    output.D1_bc=output.y_diff_unbounded_bc.D1;
    output.D2_bc=output.y_diff_unbounded_bc.D2;
    output.D4_bc=output.y_diff_unbounded_bc.D4;
    
    output.y_bc=output.y_diff_unbounded_bc.y_list;
    syms x 
    U_sym=tanh(x);
    dU_sym=diff(U_sym);
    ddU_sym=diff(dU_sym);

    output.U_bar=double(subs(U_sym,output.y_bc));
    output.d_U_bar=double(subs(dU_sym,output.y_bc));
    output.dd_U_bar=double(subs(ddU_sym,output.y_bc));      

elseif strcmp(flag.mean,'bickley_jet')    
    %conduct the coordinate transform suggested by 
    %Deloncle A, Chomaz JM, Billant P. Three-dimensional stability of a horizontally sheared flow in a stably stratified fluid. Journal of Fluid Mechanics. 2007 Jan;570:297-305.
    [y_diff_unbounded_full,output.y_diff_unbounded_bc]=channel2unbounded(y_diff_full);
    output.D1_bc=output.y_diff_unbounded_bc.D1;
    output.D2_bc=output.y_diff_unbounded_bc.D2;
    output.D4_bc=output.y_diff_unbounded_bc.D4;
    output.y_bc=output.y_diff_unbounded_bc.y_list;
    syms x 
    U_sym=sech(x)^2;
    dU_sym=diff(U_sym);
    ddU_sym=diff(dU_sym);

    output.U_bar=double(subs(U_sym,output.y_bc));
    output.d_U_bar=double(subs(dU_sym,output.y_bc));
    output.dd_U_bar=double(subs(ddU_sym,output.y_bc));         
else
    error('Wrong flag.mean');
end

for re_ind=1:length(flag.Re_list)
    Re=flag.Re_list(re_ind);
    for rotation_ind=1:length(flag.rotation_list)
        rotation=flag.rotation_list(rotation_ind);
        for Ri_ind=1:length(flag.Ri_bulk_list)
            Ri_bulk=flag.Ri_bulk_list(Ri_ind);
            if strcmp(flag.post,'eig_stratified_spanwise_tanh_asymptotic') || ...
                    (strcmp(flag.post,'eig_stratified_spanwise_bickley_asymptotic'))
                Fh=1/sqrt(Ri_bulk);
                flag.kz_list=flag.Fh_kz_list/Fh;
            end
            
            for Pr_ind=1:length(flag.Pr_list)
                Pr=flag.Pr_list(Pr_ind);
                for kx_ind=1:length(flag.kx_list)
                    kx=flag.kx_list(kx_ind);
                    for kz_ind=1:length(flag.kz_list)
                        kz=flag.kz_list(kz_ind);
                        K2= kx^2+kz^2; % Total wave number, in convienent for calculation
                        zi=sqrt(-1);
                        
                        
                        A11=(output.D4_bc-2*K2*output.D2_bc+K2^2*output.I_bc)/Re+diag(output.dd_U_bar)*zi*kx*output.I_bc-zi*kx*diag(output.U_bar)*(output.D2_bc-K2*output.I_bc); %%Orr-Sommerfeld operator
                       
                        A21= -zi*kz*diag(output.d_U_bar)*output.I_bc; %Coulping operator
                        A22= -zi*kx*diag(output.U_bar)*output.I_bc+1/Re*(output.D2_bc-K2*output.I_bc); %Squire operator  
                        %%Inverse of [Laplacian, zeros; zeros, I]
                        inv_lap=inv([output.D2_bc-K2*output.I_bc, zero_bc; zero_bc, output.I_bc]);
                        A= inv_lap*[A11, zero_bc-rotation*1i*kz*output.I_bc; A21+rotation*1i*kz*output.I_bc, A22];
                        
                        
                        dz_rho_mean=-output.I_bc; %%For spanwise stratification, I can only work on this version...
                        
                        A_stratified=[A, (-Ri_bulk)*inv_lap*[ -zi*kz*output.D1_bc; -zi*kx*output.I_bc];
                                    -dz_rho_mean*[zi*kz*output.D1_bc, zi*kx*output.I_bc]/K2,-zi*kx*diag(output.U_bar)*output.I_bc+1/Re/Pr*(output.D2_bc-K2*output.I_bc)];
                        A=A_stratified;
                        [V,D]=eig(A);
                        if isfield(flag,'store_all_eig_vec') && flag.store_all_eig_vec
                            output.eig_vec_all{re_ind,rotation_ind,Ri_ind,Pr_ind,kx_ind,kz_ind}=V;
                        end
                        output.eig_val{re_ind,rotation_ind,Ri_ind,Pr_ind,kx_ind,kz_ind}=diag(D);
                     
                        [output.eig_real_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind),ind]=max(real(diag(D)));                          
                        %output.eig_imag_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=imag(D(ind,ind));
                        % 
                        if max(real(diag(D)))>0
                            output.eig_imag_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=imag(D(ind,ind));
                        else 
                            output.eig_imag_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=NaN;
                        end
                        output.eig_vec{re_ind,rotation_ind,Ri_ind,Pr_ind,kx_ind,kz_ind}=V(:,ind);
                    end
                end
                output.eig_real_max_kx_kz(re_ind,rotation_ind,Ri_ind,Pr_ind)=...
                    max(max(output.eig_real_max{re_ind,rotation_ind,Ri_ind,Pr_ind}));
            end
        end
    end
end

% Preserve solver state needed by the existing postprocessing.
output.kz_list=flag.kz_list;
output.Re=Re;
output.kx=kx;
output.kz=kz;
end

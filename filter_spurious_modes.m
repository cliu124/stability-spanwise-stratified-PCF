function filtered=filter_spurious_modes(output,output_high_Ny,tolerance)
% Retain lower-grid eigenpairs that converge against the higher-grid spectrum.
% Matches are one-to-one, prioritizing the smallest complex relative error.
% Retained spectra are sorted by decreasing real part, then imaginary part.
if nargin<3
    tolerance=0.01;
end
assert(isequal(size(output.eig_val),size(output_high_Ny.eig_val)), ...
    'The two datasets must contain the same parameter combinations.');
assert(isfield(output,'eig_vec_all'), ...
    'Compute the lower-grid dataset with flag.store_all_eig_vec=true.');

filtered=output;
filtered.convergence_tolerance=tolerance;
filtered.eig_relative_error=cell(size(output.eig_val));
filtered.eig_match_index=cell(size(output.eig_val));
for spectrum_ind=1:numel(output.eig_val)
    low=output.eig_val{spectrum_ind}(:);
    high=output_high_Ny.eig_val{spectrum_ind}(:);
    % Implicit expansion compares every lower-grid mode with every higher-grid mode.
    relative_error=abs(low-high.')./max(abs(high.'),eps);
    relative_error(~isfinite(relative_error))=Inf;
    [candidate_low,candidate_high]=find(relative_error<tolerance);
    candidate_error=relative_error(sub2ind(size(relative_error),candidate_low,candidate_high));
    [~,order]=sort(candidate_error);
    used_low=false(numel(low),1);
    used_high=false(numel(high),1);
    match=zeros(numel(low),1);
    errors=nan(numel(low),1);
    for candidate_ind=order(:).'
        i=candidate_low(candidate_ind);
        j=candidate_high(candidate_ind);
        if ~used_low(i) && ~used_high(j)
            used_low(i)=true;
            used_high(j)=true;
            match(i)=j;
            errors(i)=candidate_error(candidate_ind);
        end
    end
    keep=find(used_low);
    [~,order]=sortrows([-real(low(keep)),-imag(low(keep))],[1,2]);
    keep=keep(order);
    filtered.eig_val{spectrum_ind}=low(keep);
    filtered.eig_vec_all{spectrum_ind}=output.eig_vec_all{spectrum_ind}(:,keep);
    filtered.eig_relative_error{spectrum_ind}=errors(keep);
    filtered.eig_match_index{spectrum_ind}=match(keep);

    [re_ind,rotation_ind,Ri_ind,Pr_ind,kx_ind,kz_ind]= ...
        ind2sub(size(output.eig_val),spectrum_ind);
    if isempty(keep)
        growth=NaN;
        frequency=NaN;
        filtered.eig_vec{spectrum_ind}=[];
    else
        growth=real(low(keep(1)));
        % Preserve the solver convention: stable modes have no plotted frequency.
        frequency=NaN;
        if growth>0
            frequency=imag(low(keep(1)));
        end
        filtered.eig_vec{spectrum_ind}=output.eig_vec_all{spectrum_ind}(:,keep(1));
    end
    filtered.eig_real_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=growth;
    filtered.eig_imag_max{re_ind,rotation_ind,Ri_ind,Pr_ind}(kx_ind,kz_ind)=frequency;
end

for parameter_ind=1:numel(filtered.eig_real_max)
    growth=filtered.eig_real_max{parameter_ind};
    growth=growth(~isnan(growth));
    if isempty(growth)
        filtered.eig_real_max_kx_kz(parameter_ind)=NaN;
    else
        filtered.eig_real_max_kx_kz(parameter_ind)=max(growth);
    end
end
end

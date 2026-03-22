This code will reproduce the results of plane Couette flow/ plane Poiseuille flow/hyperbolic tan shear flow with spanwise stratification from the paper 

Facchini G, Favier B, Le Gal P, Wang M, Le Bars M. The linear instability of the stratified plane Couette flow. Journal of Fluid Mechanics. 2018 Oct;853:205-34.

Le Gal P, Harlander U, Borcia ID, Le Dizès S, Chen J, Favier B. Instability of vertically stratified horizontal plane Poiseuille flow. Journal of Fluid Mechanics. 2021 Jan;907:R1.

Deloncle A, Chomaz JM, Billant P. Three-dimensional stability of a horizontally sheared flow in a stably stratified fluid. Journal of Fluid Mechanics. 2007 Jan;570:297-305.

Run the code main_eig_stratification_z_validation.m within this repository and change different options of 
flag.post as one of below:

% 'eig_stratified_FFLWL_figure5a_validation',...
% 'eig_stratified_FFLWL_figure10_validation',...
% 'eig_stratified_spanwise_LHBLCF_figure6b_validation',...
% 'eig_stratified_spanwise_tanh_DCB_figure4a'

Then it will reproduce the eigenvalue results in the above papers, with sample figure within this repository and also file "stability-spanwise-stratified-PCF.pdf". Note that for some results, the grid resolution in y, kx, or kz might need to be refined. 

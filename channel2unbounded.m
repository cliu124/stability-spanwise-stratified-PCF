function [y_diff_unbounded_full,y_diff_unbounded_bc]=channel2unbounded(y_diff_full)

%%This script transform the Chebyshev domain [-1,1] to the [-inf, int]
%%This is used for the unbounded background shear flow...

%%Input: 
%%y_diff_unbounded_full={y_list_full,D1_full,D2_full,D3_full,D4_full,Iw_root_full,weight_full,residue_D4_bc,residue_y_list};
%%This is the output of the grid_diff.m
%%in there, we should have the grid points for the channel, the first to
%%forth order differential matrix where no boundary condition are
%%implemented. Also, the weight_full will be used to generate the
%%integration weight for the boundary layer.

%%Output:
%%y_diff_unbounded_full={y_list_unbounded_list,D1_unbounded_full,D2_unbounded_full,D3_unbounded_full,D4_unbounded_full,Iw_root_full,weight_unbounded_full,0,0};
%%(The all derivative differential matrix do not consider boundary condition, namely they are full matrix)
%%y_list_unbounded_full: the grid points in the wall normal location for the
%%unbounded domain.
%%D1_unbounded_full: the first order derivative of the unbounded domain.
%%D2_unbounded_full: the second order derivative of the unbounded domain.
%%D3_unbounded_full: the third order derivative of the unbounded domain.
%%D4_unbounded_full: the fourth order derivative of the unbounded domain.
%%Iw_root_full: The square root of the integration weight, in the diagonal
%%place of a matrix
%%weight_unbounded_full: the first order integration weight for the unbounded domain.

%%This script relies on chebdif() of Wiediemann (2000) and Dmat() that from
%%appendix of Schmid & Henningson (2001).  Dmat() can be replaced through using the integration weight of 
%%clencurt quadrature from the weiside: https://people.maths.ox.ac.uk/trefethen/clencurt.m


%Author: Chang Liu
%Date: 2021/07/23

%%Set up default value if the user do not give.
if nargin<1 || isempty(y_diff_full)
  error('Please provide structure generated from grid_diff_cheb.');
end


%%Perform the coordinate transform.
syms sym_eta sym_x


%%x\in [-1, 1], \eta\in [-\infty,\infty]

%%The transformation between these domains are:
%%eta=x/sqrt(1+x^2) and x=eta/sqrt(1-eta^2)
%%This is suggested on p. 487 of Schmid & Henningson (2001) and the same
%%transformation is used in the following paper
%%Deloncle A, Chomaz JM, Billant P. Three-dimensional stability of a horizontally sheared flow in a stably stratified fluid. Journal of Fluid Mechanics. 2007 Jan;570:297-305.

% sym_x=(sym_eta*b-a)/(a+sym_eta);
% sym_x=sym_eta/sqrt(1-sym_eta^2);
sym_x=sym_eta/sqrt(1+sym_eta^2);
dx_deta_sym=diff(sym_x);
d2x_deta2_sym=diff(sym_x,2);
d3x_deta3_sym=diff(sym_x,3);
d4x_deta4_sym=diff(sym_x,4);

%%Below is to modify the integration weight that is used for the
%%integration.
clear sym_x sym_eta;
syms sym_eta sym_x;
sym_eta=sym_x/sqrt(1-sym_x^2);
%sym_eta=sym_x/sqrt(1+sym_x^2);

x=y_diff_full.y_list;
D1=y_diff_full.D1;
D2=y_diff_full.D2;
D3=y_diff_full.D3;
D4=y_diff_full.D4;
% eta=x./sqrt(1+x.^2);
eta=x./sqrt(1-x.^2);

dx_deta=double(subs(dx_deta_sym,eta));
d2x_deta2=double(subs(d2x_deta2_sym,eta));
d3x_deta3=double(subs(d3x_deta3_sym,eta));
d4x_deta4=double(subs(d4x_deta4_sym,eta));

%%modify the differential matrix for the boundary layer.
D1_unbounded_full=diag(dx_deta)*D1;
D2_unbounded_full=diag(dx_deta.^2)*D2+diag(d2x_deta2)*D1;
D3_unbounded_full=diag(dx_deta.^3)*D3+3*diag(dx_deta.*d2x_deta2)*D2+diag(d3x_deta3)*D1;
D4_unbounded_full=diag(dx_deta.^4)*D4...
    +6*diag(dx_deta.^2.*d2x_deta2)*D3...
    +4*diag(d3x_deta3.*dx_deta)*D2...
    +3*diag(d2x_deta2.^2)*D2+diag(d4x_deta4)*D1;


%%Update 2020/03/05, This D4_BL_full needs to be updated. This should
%%follow Chapter 14 of Trefethen (2000) Spectral method in MATLAB. So the
%%D4_BL_full operator will directly involve the boundary condition of the
%%first order derivative. Previous eigenvalue computation and so on are not
%%using this. However, if we would like to compute the H2 norm of the
%%boundary layer (wall normal velocity and vorticity form), 
%it is necessary to implement the boundary conditions of this correct...
%This will overwrite the old D4_BL_full. The default ymin=0
%p(eta)=(eta-eta_min)*(eta_max-eta)*q(eta)
%d4p=-4*(eta-eta_min)*d3q-4*(eta-eta_max)d3q-12*d2q-(eta-eta_min)*(eta-eta_max)d4q

%%%This fourth derivative boundary condition is assumed to be satisfied...
%%This is because following the traditional boundary condition
%%implementation... it will only left the fourth order term and other terms
%%are become zero....

% ymin=-Inf; ymax=Inf;
% S=diag([0;1./(eta(2:end-1)-ymin)./(ymax-eta(2:end-1));0]);
% D4_unbounded_full=(-4*diag(eta-ymin)*D3_unbounded_full-4*diag(eta-ymax)*D3_unbounded_full...
%             -12*D2_unbounded_full-diag((eta-ymin).*(eta-ymax))*D4_unbounded_full)*S;

y_list_unbounded_list=eta;

%%Set up the output
y_diff_unbounded_full.y_list=y_list_unbounded_list;
y_diff_unbounded_full.D1=D1_unbounded_full;
y_diff_unbounded_full.D2=D2_unbounded_full;
y_diff_unbounded_full.D3=D3_unbounded_full;
y_diff_unbounded_full.D4=D4_unbounded_full;

y_diff_unbounded_bc.y_list=y_list_unbounded_list(2:end-1);
y_diff_unbounded_bc.D1=D1_unbounded_full(2:end-1,2:end-1);
y_diff_unbounded_bc.D2=D2_unbounded_full(2:end-1,2:end-1);
y_diff_unbounded_bc.D3=D3_unbounded_full(2:end-1,2:end-1);
y_diff_unbounded_bc.D4=D4_unbounded_full(2:end-1,2:end-1);


    
end



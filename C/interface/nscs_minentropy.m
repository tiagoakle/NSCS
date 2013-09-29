%This script calls the nscs solver to solve an entropy minimization problem

[ret,warn] = loadlibrary('../C/Lib/nscs','../C/Include/matlab_include/matlab_lib.h');
% ENTROPY PROBLEM

clear; clc;
 
M = 20;
N = 50;
 
A    = randn(M,N);
xx   = 3*ones(N,1);
b    = A*xx;
d    = ones(N,1);
 
pars.echo   = 4;
pars.beta   = 0.99;
pars.trace  = 3;
pars.secord = 0;

AA = [[sparse(N,N),speye(N),sparse(N,N)];[sparse(M,2*N),A]];
bb = [ones(N,1);b];
c  = [-ones(N,1);zeros(2*N,1)];

%Define the permutation to put A into the ordering 
% which is appropriate for nscs
permute = zeros(3*N,1);
permute(1:3:3*N) = [1:N];
permute(2:3:3*N) = N+[1:N];
permute(3:3:3*N) = 2*N+[1:N]; 

tK = 3*ones(N,1);
nK = 3*ones(N,1);
k_count = N;


AA = AA(:,permute);
[AI,AJ,AV] = find(AA);
nnzA = size(AI,1);
m    = M+N;
n    = 3*N;
c    = c(permute);

% starting point:
u0  = -ones(N,1);
v00 = ones(N,1);  
x0  = 0.5*ones(N,1);
xx0  = [u0;v00;x0];
xx0  = xx0(permute);

%Allocate the memory for the work vectors
p_y  = libpointer('doublePtr',zeros(m,1));
p_x  = libpointer('doublePtr',xx0);
p_t  = libpointer('doublePtr',0);
p_s  = libpointer('doublePtr',zeros(n,1));
p_k  = libpointer('doublePtr',0);
wy   = 0;
ws   = 0;
wk   = 0;
wt   = 0;

% stopping constants:
relstopP  = max(1,norm([AA,bb],'inf'));
relstopD  = max(1,norm([AA',speye(n),-c],'inf'));
relstopG  = max(1,norm([-c',bb',1],'inf'));

%Parameters from the default of coneopt
max_iter  = 500;
max_center_iter  = 50;
theta     = 0.8;
lscaff    = 0.94;
lsccent   = 0.5;
eta       = 0.9995;
max_backtrack = 300;
delta     = 1.e-10; 
gamma     = 1.e-10;
beta      = 0.99;
p_relstop =  max(1,norm([AA,bb],'inf'));
d_relstop = max(1,norm([AA',speye(n),-c],'inf'));
g_relstop =  max(1,norm([-c',bb',1],'inf'));
p_rho = 1e-6; 
d_rho = 1e-6;
a_rho = 1e-6; 
rho_g = 1e-6;
rhoI  = 1e-6;
rhoM  = 1e-6;

ret = calllib('nscs','nscs_no_structs',k_count,nK,tK,AI,AJ,AV,nnzA,m,n,...
                      c,b,...
                      p_y,p_x,p_s,p_t,p_k,...
                      wy,ws,wt,wk,...
                      max_iter,max_center_iter,  theta,  lscaff,  lsccent,...
                      eta, max_backtrack,  delta,  gamma,  beta,...
                      p_relstop, d_relstop, g_relstop,...
                      p_rho, d_rho, a_rho, rho_g,...
                      rhoI,  rhoM);



%void nscs_no_structs(csi k_count, csi* nK, int* tK,\
%                     csi *AI, csi* AJ, double* AV, csi nnzA, csi m, csi n,\
%                     double *c, double *b,\
%                     double* y, double* x, double* s, double* tau, double* kappa,\
%                     int wy, int ws, int wt, int wk,\
%                     int max_iter, int max_center_iter, double theta, double lscaff, double lsccent,\
%                     double eta, int max_backtrack, double delta, double gamma, double beta,\
%                     double p_relstop, double d_relstop, double rel_gap_relstop,\
%                     double p_rho, double d_rho, double a_rho, double rho_g,\
%                     double rhoI, double rhoM);


unloadlibrary nscs;

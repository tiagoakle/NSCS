clear all
%Runs the GP instance read by read_gp.m in coneopt
addpath '../../';
clear all

%Read a gp
file_name = './gp/beck751.eo';
[AA,bb,cc,num_ter,num_var,num_con] = read_gp(file_name);

% starting point:
t00 = 1;
up0 = ones(num_var,1);
um0 = ones(num_var,1);
w00 = -ones(num_ter,1);
v00 = ones(num_ter,1);  
y00 = 0.5*ones(num_ter,1);

%Extract the problem data and build the problem structure
problem = struct;
problem.A = AA;
problem.b = bb;
problem.c = cc;

%Problem parameters
problem.m =2*num_ter + num_con+1;
problem.n = 3*num_ter + 2*num_var + 1;
problem.n_free = 0;
problem.n_constrained = 2*num_var+1+3*num_ter;
problem.n_pos       = 2*num_var+1;
problem.soc_cones   = 0;
problem.n_soc_cones = 0;
problem.n_sdp_cones = 0;
problem.sdp_cones     = 0;
problem.n_exp_cones   = num_ter;
problem.n_power_cones = 0;

x0c  = [t00;up0;um0;w00;v00;y00];
x0f  = [];
nscs 


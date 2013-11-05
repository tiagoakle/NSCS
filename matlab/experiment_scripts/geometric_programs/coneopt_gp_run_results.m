clear all
%Runs the GP instance read by read_gp.m in coneopt
addpath '../../../coneopt';
clear all
read_gp

% build cone:
K.npos = 2*num_var+1;
K.npow = 0;
K.nexp = num_ter;
K.nlog = 0;
K      = getbarrpar(K);

%Set the parameters
pars.n = 3*num_ter + 2*num_var + 1;
pars.m = 2*num_ter + num_con+1;
pars.echo = 4;
pars.secord = 1;
pars.centmeastype = 5;

% starting point:
t00 = 1;
up0 = ones(num_var,1);
um0 = ones(num_var,1);
w00 = -ones(num_ter,1);
v00 = ones(num_ter,1);  
y00 = 0.5*ones(num_ter,1);

v0.x  = [t00;up0;um0;w00;v00;y00];

% call to coneopt:
R = coneopt(AA,bb,cc,v0,K,pars);


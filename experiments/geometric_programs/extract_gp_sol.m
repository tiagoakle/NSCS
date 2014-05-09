%Given a solution to the GP calculated by NSCS extract the original variables and evaluate
%the feasibility and objective
function [x_eo,obj,feas] = extract_gp_sol(y,xc,s,t,k,AA,bb,cc,num_ter,num_var,num_con)
%function [AA,bb,cc,num_ter,num_var,num_con] = read_gp(file_name) 
%Extract the variables 
x = xc/t;
x_eo = x(1+1:1+num_var)-x(1+num_var+1:2*num_var+1);
%Evaluate the objective
u    = x(1+2*num_var+1:1+2*num_var+num_ter)
v    = x(1+2*num_var+1+num_ter:1+2*num_var+2*num_ter)
w    = x(1+2*num_var+1+2*num_ter:1+2*num_var+3*num_ter)

%Feasibility stuff
w.*exp(u./w)-v

end

%Reads a GP in the format provided by Erling
%And sets it up in the form solved by coneopt
% [0,A,-A,-I, ,  ] 
% [0,0,0 ,0 ,0, I] 
% [1,0,0 ,0 ,ix_0] 
% [0,0,0 ,0 ,IX  ] 
%[obj_val, x_plus,x_minus,u,v,w]

%Each "term" gives rise to an exponential cone and a triplet u,v,w

file_name = './gp/beck751.eo';
f         = fopen(file_name,'r');
%Read the number of constraints
num_con   = fscanf(f,'%i',1);
%Read the number of variables
num_var   = fscanf(f,'%i',1);
%Read number of terms 
num_ter   = fscanf(f,'%i',1);
%Read the 'c' coefficients for the terms
[c_coef,c_coef_count] = fscanf(f,'%g',num_ter);
%Do some checking
if(c_coef_count ~= num_ter); error('Unable to read c coefficients'); end
%Read the constraint indices
[constraint_index,constraint_index_count] = fscanf(f,'%i',num_ter);
if(constraint_index_count~=num_ter); error('Unable to read the constraint indices'); end
%Not an efficient implementation
I = [];
J = [];
V = []; 

%Define the vectors that indicate which terms belong to each constraint
constraints = sparse(constraint_index+1,[1:num_ter]',ones(num_ter,1),num_con+1,num_ter,num_ter);

%Define a vector for the objective
a_0 = zeros(3*num_ter+num_var,1);
while true
    row = fscanf(f,'%i %i %e',3);
    if isempty(row)
        break;
    end
    %The data is 0 indexed
    I = [I;row(1)+1];
    J = [J;row(2)+1];
    V = [V;row(3)];
end
neg_1 = sparse(num_con+1,1);
neg_1(1) = -1;
neg_0 = sparse(num_ter,1); %Column of zeros
A   = sparse(I,J,V,num_ter,num_var);

AA  = [[neg_0,A,-A,-speye(num_ter),sparse(num_ter,2*num_ter)];
      [neg_1,sparse(num_con+1,2*num_var+num_ter),constraints,sparse(num_con+1,num_ter)];
      [sparse(num_ter,2*num_var+2*num_ter+1),speye(num_ter)]];

bb  = zeros(2*num_ter+num_con+1,1);
bb(1:num_ter) = -log(c_coef);
bb(num_ter+1) = 0;
bb(num_ter+2:2*num_ter+num_con+1) = 1;

cc = zeros(3*num_ter+2*num_var+1,1);
cc(1) = 1;

fclose(f);

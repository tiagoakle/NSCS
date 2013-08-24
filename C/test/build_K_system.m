%Reads the A and H matrices and b,c vectors from the minentropy system
% forms the KKT system 
%      [ dI     A]
% K =  [A' -mH-gI], with d=delta=1.e-10 g=gamma=1.e-10 mu=1
% Forms rhs = [b;-c];
% Solves K*sol=rhs

%Saves K in matlab_minentropy_K, sol in matlab_minentropy_sol.csv
%Saves A in minentropy_A.csv and H in minentropy_H.csv
%Saves rhs in minentropy_rhs.csv
%Saves b in minentropy_b.csv
%Saves c in minentropy_c.csv
%Saves r1,r2,r5 in minentropy_r1,...,minentropy_r5
%Saves [mu,tau,kappa,r3,r4] in minentropy_doubles;

load './test_data/minentropy_system.mat'
[m,n] = size(A);

%Choose the regularization paramters and centrality value
mu = 1;
gamma = 1.e-10;
delta = 1.e-10;

nnzA       = nnz(A);
nnzH       = nnz(H);
nnzK       = m + 2*nnzA+nnzH;

K  = [[delta*speye(m,m),A];[A',-mu*H-gamma*speye(n,n)]];
rhs= [b;-c];

%Solve the system
sol = K\rhs;

doubles =[mu,tau,kappa,r3,r4]'; 
%Save the matrices and vectors
write_vector_to_csv('./test_data/minentropy_sol.csv',sol);
write_vector_to_csv('./test_data/minentropy_rhs.csv',rhs);
write_vector_to_csv('./test_data/minentropy_b.csv',b);
write_vector_to_csv('./test_data/minentropy_c.csv',c);

write_vector_to_csv('./test_data/minentropy_r1.csv',r1);
write_vector_to_csv('./test_data/minentropy_r2.csv',r2);
write_vector_to_csv('./test_data/minentropy_r5.csv',r5);
write_vector_to_csv('./test_data/minentropy_doubles.csv',doubles)

write_matrix_to_csv('./test_data/matlab_built_K.csv',K);
write_matrix_to_csv('./test_data/minentropy_A.csv',A);
write_matrix_to_csv('./test_data/minentropy_H.csv',H);

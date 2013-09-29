%Builds a matrix A, a vector c and a vector b 
%to be used in a nscs_minentropy example

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

% stopping constants:
relstopP  = max(1,norm([AA,bb],'inf'));
relstopD  = max(1,norm([AA',speye(n),-c],'inf'));
relstopG  = max(1,norm([-c',bb',1],'inf'));

%Write A to test_data/minentropy/A.bin
write_matrix_bin('./test_data/minentropy/A.bin',AA);
%Write b and c to test_data/minenetropy/b.bin and c.bin
write_vector_bin('./test_data/minentropy/c.bin',c);
write_vector_bin('./test_data/minentropy/b.bin',b);
%Write vector with relstopP relstopD relstopG
params = [relstopP,relstopD,relstopG];
write_vector_bin('./test_data/minentropy/relstop.bin',params);
write_vector_bin('./test_data/minentropy/x0.bin',xx0);

save('./test_data/minentropy/matlab_copy');

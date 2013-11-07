%Runs the list of problems in standard form with the 3 solvers 
%The regular mehrotra lp code
%The non_symmetric_long_step with second order
%The non_symmetric_long_step with second order and arc_search

clear all
  
addpath '../matlab'
addpath '../coneopt'
%Load the file that contains the indices for the
%ufget netlib lps which are in standard form
load 'standard_form_indices.mat' 

problem_count = length(st_ix);

fid = 1;
fprintf(fid,'Will benchmark against %i problems\n',length(st_ix));

results = {{'PNAME','MLS','MLSSO','MLSAS','m','n'}};
for j = 1:problem_count
    %Get the ufget id 
    problem_uf_ix = st_ix(j);
    %Get the problem from ufget
    P = UFget(problem_uf_ix);  
    %Extract the name
    prob_name = [P.name];
    %Substitute front slash for space
    prob_name(find(prob_name=='/'))=' ';
    
    %Extract the problem data and build the problem structure
    A = P.A;
    b = P.b;
    c = P.aux.c;
    [m,n] = size(A);
    %Call the solvers 
    [x,y,s,info] = mehrotra_lp_solver(A,b,c);
    [x_ns,y_ns,s_ns,info_ns_se] = non_symmetric_long_step(A,b,c,'secord');
    [x_ns,y_ns,s_ns,info_ns_as] = non_symmetric_long_step(A,b,c,'arc_search');
    
    %Append the new results 
    results = {results{:},{prob_name,info.iter,info_ns_as.iter,info_ns_se.iter,m,n}};

    %Print the results so far
    fprintf('%25s %6s %6s %6s %4s %4s\n',results{1}{1},results{1}{2},results{1}{3},...
                                results{1}{4},results{1}{5},results{1}{6});
    for k = 2:j+1
        fprintf('%25s   %3i    %3i    %3i  %4i  %4i\n',results{k}{1},results{k}{2},...
                                         results{k}{3},results{k}{4},...
                                         results{k}{5},results{k}{6});
    end
    %Save the present progress
    save('mehrotra_lpnetlib_benchmark.mat','results');
end

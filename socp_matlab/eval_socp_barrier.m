%Evaluates the gradient of the socp barrier for a product of cones
function f = eval_socp_barrier(problem,xc) 
    f = 0; 
    i = 1; %Index of the first variable of the cone
    for(j=1:problem.n_soc_cones)
       %Calculate the barrier for this cone
       f = -log(xc(i)^2-norm(xc(i+1:i+problem.soc_cones(j)-1))^2)+f;
       i = i+problem.soc_cones(j);
    end 
end



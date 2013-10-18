%Evaluates the gradient of the socp barrier for a product of cones
function g = eval_socp_gradient(problem,xc) 

    g = zeros(size(xc));
    i = 1; %Index of the first variable of the cone
    for(j=1:problem.n_soc_cones)
       %Calculate 2/x'Jx
       xjxi = xc(i)^2-norm(xc(i+1:i+problem.soc_cones(j)-1))^2;
       xjxi = -1/xjxi;
       g(i:i+problem.soc_cones(j)-1) = xjxi*[xc(i);-xc(i+1:i+problem.soc_cones(j)-1)];
       i = i+problem.soc_cones(j);
    end 
end



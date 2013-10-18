%Evaluates the hessian of the socp barrier for a product of cones
function H = eval_socp_hessian(problem,xc) 
    nnz = sum(problem.soc_cones.^2);
    siz = sum(problem.soc_cones);
    H   = sparse([],[],[],siz,siz,nnz); 
    
    i = 1; %Index of the first variable of the cone
    for(j=1:problem.n_soc_cones)
       %Calculate 1/x'Jx
       xjxi = xc(i)^2-norm(xc(i+1:i+problem.soc_cones(j)-1))^2;
       xjxi = 1/xjxi;
       Jx   = [xc(i);-xc(i+1:i+problem.soc_cones(j)-1)];
       H_s  = 2*xjxi^2*Jx*Jx'-xjxi*diag([1;-ones(problem.soc_cones(j)-1,1)]);
       H(i:i+problem.soc_cones(j)-1,i:i+problem.soc_cones(j)-1) = H_s;
       i = i+problem.soc_cones(j);
    end 
end



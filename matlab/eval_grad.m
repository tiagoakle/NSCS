function [g] = eval_grad(problem,xc)
             g = zeros(problem.n_constrained,1);
             %Evaluate the positive orthant variables
             g(1:problem.n_pos) = -1./xc(1:problem.n_pos);

             %Evaluate the soc_cone gradient
             ix = problem.n_pos+1;
             for(i=1:problem.n_soc_cones)
                nrm = norm(xc(ix+1:ix+problem.soc_cones(i)-1));
                scaling = -2/(-nrm^2+xc(ix));
                g(ix)   = scaling*xc(ix);
                g(ix+1:ix+problem.soc_cones(i)-1) = -scaling*xc(ix+1:ix+problem.soc_cones(i)-1);
                ix = ix + problem.soc_cones(i);
             end

             %Evaluate the sdp_cone_gradient
             %These are the inverses of the matrices
             for(i=1:problem.n_sdp_cones)
                 M  = reshape(xc(ix:ix+sdp_cones(i)^2-1),sdp_cones(i),sdp_cones(i));
                [R] = chol(M);
                R_inv = R\eye(sdp_cones(i));
                M   = R_inv*R_inv';
                g(ix:ix+sdp_cones(i)^2-1) = M(:);
                ix = ix + problem.sdp_cones(i)^2;
             end
            
             %Evaluate the exponential cones 
             %in the pattern [x1s, x2s, x3s]
             if(problem.n_exp_cones>0)
                x1       = xc(ix:ix+problem.n_exp_cones-1);
                x2       = xc(ix+problem.n_exp_cones:ix+2*problem.n_exp_cones-1);
                x3       = xc(ix+2*problem.n_exp_cones:ix+3*problem.n_exp_cones-1);

                logx2 = log(x2);
                logx3 = log(x3);
                x2m1  = 1./x2;
                x3m1  = 1./x3;
                tmp1  = logx2-logx3;
                psi   = x3.*tmp1 - x1;
                psim1 = 1./psi;
                xi    = tmp1 - 1;


                 tmp2  = [ psim1;...
                     -x2m1.*(x3.*psim1 + 1);...
                     -xi.*psim1 - x3m1];
                g(ix:ix+3*problem.n_exp_cones-1) = tmp2(:);
            end

end

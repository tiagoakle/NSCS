%Feb 18 2013
%Plots the potential reduction along the present direction
%Requires x,y,z and d_x d_y d_z, alph_xz 
%and rho to be defined in the workspace
function [phi,phi_mod,phi_frac] = plot_present_point_and_direction(x,y,z,x_d,y_d,z_d,samples,A,b,c,max_alph_xz,rho,local_descent,eta,merit_a) 
    
    %Allocate space for the samples
    phi           = zeros(samples,1);
    phi_mod       = zeros(samples,1);
    linear        = zeros(samples,1);
    lin_pred      = zeros(samples,1); 
    %Calculate the residuals
   
    p_r  = b - A*x;
    d_r  = c - A'*y - z;
    n_p_r = norm(p_r);
    n_d_r = norm(d_r);
    gap   = abs(c'*x-b'*y);
    

    %Calculate the gradient of the objective   
    l1    = 0.5*n_p_r^2 + 0.5*n_d_r^2 + x'*z;
    g_merit = [rho/l1*(-A'*p_r+z)-1./x;rho/l1*(-A*d_r);rho/l1*(-d_r + x)-1./z];

    l2    = sum(log(x)) + sum(log(z));
    merit = rho*log(l1) - l2;
    
    local_descent = g_merit'* [x_d;y_d;z_d];
    %--------------------------------------------------------------------------------------
    %Sample the potential reduction functions
    %--------------------------------------------------------------------------------------

   %Save the present step    
    xs = x;
    ys = y;
    zs = z;
        
   fprintf('Will sample the potential functions along the search direction \n'); 
   fprintf('Local Descent prediction %d\n',local_descent);

    ix = 1;

    r= [c-z-A'*y;
        b-A*x];
    Eds = [A'*y_d+z_d;A*x_d];

    cr    = 0.5*norm(r)^2;
    br    = -r'*Eds;
    ar    = 0.5*norm(Eds)^2;
 
    for al = 0:max_alph_xz/samples:max_alph_xz-max_alph_xz/samples
        x  = xs+al*x_d;
        y  = ys+al*y_d;
        z  = zs+al*z_d;
        %calculate the potential function
        p_r= b-A*x;
        d_r= c-z-A'*y;
        phi(ix)     = rho*log(al^2*ar+al*br+cr+x'*z)-sum(log(x.*z));
        phi_mod(ix) = rho*log(0.5*norm(p_r)^2+0.5*norm(d_r)^2+x'*z) - sum(log(x)) - sum(log(z));
        linear(ix)   = merit_a + eta*local_descent*al;
        lin_pred(ix) = merit_a + local_descent*al;
        ix = ix + 1;
   end
    fprintf('Sampled the potential reduction functions\n');

%figure
%plot(phi);
%title('Potential reduction ye');

figure
hold on;
plot(phi);
plot(phi_mod,'k');
plot(linear,'g');
plot(lin_pred,'r');
title('Potential reduction');
legend('function poly','function','linear condition','linear prediction');
hold off;

end

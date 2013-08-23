function [d,CF] = solve_linear_system(H,mu,A,b,c,tau,kappa,r1,r2,r3,r4,r5,pars)
    
    %Debug flag 
    print_t = false;
    print_c = true;
    %Get the size of the problem
    [m,n] = size(A);
    
    [slv_aug,h_free] = c_factor();
     mat_slv_aug           = factor();

    %Eventually this will be implemented in C
    %Solves  [    A -b      ] dy   r1
    %        [-A'    c -I   ] dx   r2
    %        [b' -c'     -1 ] dt = r3
    %        [    mH   I    ] ds   r5
    %        [       k    t ] dk   r4

    %Forms 
    %        [    A -b      ] dy   r1
    %        [-A'    c -I   ] dx   r2
    %        [b' -c'     -1 ] dt = r3
    %        [    mH   I    ] ds   r5
    %        [       h   1  ] dk   r6
    %with h = t/k r6 = r4/kappa

    h  = kappa/tau;
    r6 = r4/tau;

    %Forms 
    %        [    A  -b      ] dy     r1
    %        [A' -mH -c      ] dx     r7
    %        [b' -c'  h      ] dt =   r8
    %        [   mH     I    ] ds     r5
    %        [        h   1  ] dk     r6
    % with r7 = -r2-r5; r8 = r3+r6;
    
    r7 = -(r2+r5);
    r8 = r3+r6;

    %  
    %        [    A  -b     ] dy     r1
    %        [A' -mH -c     ] dx     r7
    %        [       h_2    ] dt =   r9
    %        [   mH     I   ] ds     r5
    %        [       h    1 ] dk     r6
      
    %[b' -c'] [     A]^{-1}[-b]     =  h_1b
    %         [A' -mH]     [-c] 

    %[b' -c'] [     A]^{-1}[r1]     =  r_7b
    %         [A' -mH]     [r7] 

    %We solve 
    % tm_1 =  [     A]^{-1}[b ]    
    %         [A' -mH]     [-c] 
    % and keep it to calculate h_1b and r_7b more efficiently
    if(print_t) fprintf('Norm [b;-c] %g\n',norm([b;-c])); end;
    tm_1   = slv_aug([b;-c]); 
    if(print_t) fprintf('Norm tm1 %g\n',norm(tm_1)); end;
    if(print_c) tm1b = mat_slv_aug([b;-c]); fprintf('Difference between solutions to tm %g\n', norm(tm_1-tm1b)); end;
    h_1b   = tm_1'*[-b;-c];
    r_7b   = tm_1'*[r1;r7];
 
    r9  = r8 - r_7b;
    h_2 = h-h_1b;
    
    %We now start the back substitution
    dt  = r9/h_2;
    %reuse tm_1
    tm_1 = slv_aug([r1;r7]-dt*[-b;-c]);
    %Extract dy dx from dt
    dy   = tm_1(1:m);
    dx   = tm_1(m+1:m+n);

    %Now back substitute to get ds and dkappa
    ds   = r5-mu*H*dx;
    dk   = r6-h*dt;

    %Check the residuals
    n_res_1 = norm(A*dx-dt*b-r1);
    n_res_2 = norm(-A'*dy + dt*c -ds -r2);
    n_res_3 = norm(b'*dy-c'*dx -dk-r3);
    n_res_5 = norm(mu*H*dx+ds-r5);
    n_res_4 = norm(kappa*dt+tau*dk-r4);

    if(print_t) fprintf('Residuals r1 %g, r2 %g, r3 %g, r5 %g, r4 %g \n',n_res_1,n_res_2,n_res_3,n_res_5,n_res_4); end;
    %Solves  [    A -b      ] dy   r1
    %        [-A'    c -I   ] dx   r2
    %        [b' -c'     -1 ] dt = r3
    %        [    mH   I    ] ds   r5
    %        [       k    t ] dk   r4


    %PAck the solution into a struct
%    v.dx     = d{1};
%    v.dtau   = d{2};
%    v.dy     = d{3};
%    v.ds     = d{4};
%    v.dkappa = d{5};
  
    d    = {}; 
    d{1} = dx;
    d{2} = dt;
    d{3} = dy;
    d{4} = ds;
    d{5} = dk;
    d{6} = 0;  %This is some kind of error flag??
    
    %This variable is used to accelerate the BFGS heuristic 
    CF = [];
    
    %Free the factorization
    h_free();
  
    function [h_solver,h_free] = c_factor()
        %Allocate the variables that will hold 
        %the factorization 
        lppI     = libpointer('int32Ptr',1);
        lppP     = libpointer('int32Ptr',1);
        lppV     = libpointer('doublePtr',1.0);
        numeric  = libpointer('voidPtr');
        
        %Convert H and A to triplet form
        [aI,aJ,aV] = find(A);
        nnzA       = size(aI,1);
        [hI,hJ,hV] = find(H);
        nnzH       = size(hI,1);
        finalnnz   = m + 2*nnzA+nnzH;
        
        delta = 1.e-10;
        gamma = 1.e-10;
        
        %Do this to make sure that the 
        %representation is actually updated 
        hI(1:end) = hI(1:end)-1;
        hJ(1:end) = hJ(1:end)-1;
        aI(1:end) = aI(1:end)-1;
        aJ(1:end) = aJ(1:end)-1;

          %Call the method that builds the CCO kkt system 
          ret     = calllib('liblinear_solvers','form_kkt_system',mu,hI,hJ,hV,nnzH,aI,aJ,aV,nnzA,m,n,delta,gamma,lppI,lppP,lppV);
          fprintf('ret from form_kkt_system %i\n',ret);
          %Call the method that makes the factorization
          ret     = calllib('liblinear_solvers','factor_kkt_system',numeric,lppI,lppP,lppV,n+m);
          fprintf('ret from factor_kkt_system %i\n',ret);
          h_solver = @solve;
          h_free   = @free;

          function x = solve(rhs)
             %Call the function that solves the system 
            lpX     = libpointer('doublePtr',zeros(size(rhs,1),1));

            %Make sure the vector is actually formed
            %this matlab absurdity is nuts
            rhs2 = zeros(size(rhs));
            for iter = 1:size(rhs,1)
                rhs2(iter) = rhs(iter);
            end

            ret     = calllib('liblinear_solvers','solve_factored_system',numeric,lppI,lppP,lppV,rhs2,lpX);
            fprintf('ret from sove_factored_system %i\n',ret);
            fprintf('Norm of the solution %g\n', norm(lpX.value));
            x       = lpX.Value;
          end 

          function free()
            %Function that calls umfpack to free the variables
            ret     = calllib('liblinear_solvers','free_factorization',numeric,lppI,lppP,lppV);
          end

    end


    function h_solver = factor()
       %Factorizes
        %[dI  A    ]    = PLDL'P'
        %[A' -H -gI]
        delta = 1.e-10;
        gamma = 1.e-10;
        [L,D,P] = ldl([[delta*speye(m,m),A];[A',-mu*H-gamma*speye(n,n)]]);
        %Instantiate a handle to the solver
        h_solver = @(y)solve(L,D,P,y);
        if(print_t) fprintf('Sizes m: %i n %i \n',m,n); end
        linopts_lt = struct;
        linopts_ut = struct;
        linopts_lt.LT = true;
        linopts_ut.UT = true;

        d = diag(D);
        if(print_t) fprintf('min(d): %g, max(d): %g \n',full(min(d)),full(max(d))); end;
        if(print_t) fprintf('min diag L: %g, max |L|: %g \n', min(abs(full(diag(L)))),full(max(max(abs(L))))); end;

        function y = solve(L,D,P,y)
        
            %For debug 
            y_t = y;
            %PLDL'P'x = b 
            if(print_t) fprintf('Size rhs: %i\n', size(y,1)); end;
            if(print_t) fprintf('Norm rhs: %g\n', norm(y)); end;
            y  = (P*(L'\(D\(L\(P'*y)))));
           % d = diag(D);
           % y = P'*y;
           % y = L\y;
           % y = y./d;
           % y = (L')\y;
           % y = P*y;
            if(print_t) fprintf('Norm solution %f\n',norm(y)); end;
            
            n_res = norm([[delta*speye(m,m),A];[A',-H-gamma*speye(n,n)]]*y-y_t);
            if(print_t) fprintf('N res %g\n',n_res/norm(y_t)); end;

        end

    end
       
 end


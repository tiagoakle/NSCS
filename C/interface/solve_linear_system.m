function [d,CF] = solve_linear_system(H,mu,A,b,c,tau,kappa,r1,r2,r3,r4,r5,pars)
    %Debug flag 
    print_t = false;
    %Get the size of the problem
    [m,n] = size(A);
    
    slv_aug = factor();

%`[d,CF] = slvhomkkt(v.F{3},v.mu,pars.A,pars.b,pars.c,v.tau,v.kappa,...r1,r2,r3,r4,r5,pars);
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

    if(print_t) fprintf('Size [b;-c]: %i\n',size([b;-c],1)); end;
    tm_1   = slv_aug([b;-c]); 
    if(print_t) fprintf('Norm tm1 %f\n',norm(tm_1)); end;
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

    if(print_t) fprintf('Residuals r1 %f, r2 %f, r3 %f, r5 %f, r4 %f \n',n_res_1,n_res_2,n_res_3,n_res_5,n_res_4); end;
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
        if(print_t) fprintf('min(d): %f, max(d): %f \n',full(min(d)),full(max(d))); end;
        if(print_t) fprintf('min diag L: %f, max |L|: %f \n', min(abs(full(diag(L)))),full(max(max(abs(L))))); end;

        function y = solve(L,D,P,y)
        
            %For debug 
            y_t = y;
            %PLDL'P'x = b 
            if(print_t) fprintf('Size rhs: %i\n', size(y,1)); end;
            if(print_t) fprintf('Norm rhs: %f\n', norm(y)); end;
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


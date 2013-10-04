function R = coneopt(varargin)

[v,K,pars,R]  = init(varargin{:});

while true
    
    if R.stop, break; end                % if R.stop, break loop
    
    v.F   = BarrFuncP(v.x,K,[1,1,1]);    % evaluate barrier
    
    [v,R] = slvatd_coneopt(v,K,pars,R);  % solve for atd

    if R.stop                            % stop if rounding
        R = checkstopcrit(v,pars,R);     % blocks further progress
        break;                              
    end             
        
    [v.a,R] = linesearch(v,K,pars,R);    % line search along atd
    v       = updatevars(v,K);           % take step

    v.kappa    
    [v,R]   = getresiduals(v,pars,R);    % new residuals
    R       = checkstopcrit(v,pars,R);   % check stopping crits
    
    [v,R]   = centering(v,K,pars,R);     % centering (correction)
    v.s     = -pars.A'*v.y + ...         % pd-lifting
        pars.c*v.tau - v.rD0;               
    
    v.k     = v.k + 1;                   % increment counter
    R       = gettrace(v,pars,R);        % get current data
    
    printfunc('iter',v,K,pars,R);        % print iteration info
    
end

[v,K,pars,R]  = deinit(v,K,pars,R);
printfunc('final',v,K,pars,R);

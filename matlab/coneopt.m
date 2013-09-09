function R = coneopt(varargin)

[v,K,pars,R]  = init(varargin{:});

addpath '../C/interface'  %Add the diretory with the interface m files
load_c_library()

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
    
    [v,R]   = getresiduals(v,pars,R);    % new residuals
    R       = checkstopcrit(v,pars,R);   % check stopping crits
    
    if R.stop
        fprintf('Stop command issued\n');    
    end
%    tic 
    [v,R]   = centering(v,K,pars,R);     % centering (correction)
%    center_time = toc;
%    fprintf('Centering time %i\n',center_time);

     v.s     = -pars.A'*v.y + ...         % pd-lifting
            pars.c*v.tau - v.rD0;               

    v.k     = v.k + 1;                   % increment counter
    R       = gettrace(v,pars,R);        % get current data
    
    printfunc('iter',v,K,pars,R);        % print iteration info
    
end

[v,K,pars,R]  = deinit(v,K,pars,R);
printfunc('final',v,K,pars,R);

unload_c_library()

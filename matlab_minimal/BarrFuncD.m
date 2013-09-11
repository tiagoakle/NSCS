function F = BarrFuncD(x,K,want)

% want(1): F{1} = Function value
% want(2): F{2} = Gradient
% want(3): F{3} = Hessian
%          F{4} = Feasibility (0 for feas, -1 for infeas)

%% (DUAL) Positive cone variables

xpos = x(K.posc.idx.x);

if ~isempty(xpos)
    
    [FPOS,GPOS,HPOS,feasPOS] = BarrFuncPos(xpos,K,want);
    
    % this is specific for dual cone:
    FPOS = FPOS - K.npos;
    %see e.g. p.328 bottom of NestTodd, Siam J. Optim 1998

    if feasPOS < 0
        F{4} = -1;
        return
    end
    
else
    FPOS = [ ]; GPOS = [ ]; HPOS = [ ];
end


%% (DUAL) Power cone variables

x1  = x(K.powc.idx.x1)';
x2  = x(K.powc.idx.x2)';
x3  = x(K.powc.idx.x3)';

if ~isempty(x1)
    
    %scale vars to dual: (see p. 128-130 in Chares)
    x2  = x2 ./ K.powc.a1d(1,:);
    x3  = x3 ./ K.powc.a1d(2,:);
    
    [FPOW,GPOW,HPOW,feasPOW] = BarrFuncPow(x1,x2,x3,K,want);

    % feasibility:
    if feasPOW < 0
        F{4} = -1;
        return;
    end

else
   FPOW = [ ]; GPOW = [ ]; HPOW = [ ]; 
end



%% (DUAL) Exponential cone variables

x1  = x(K.expc.idx.x1)';
x2  = x(K.expc.idx.x2)';
x3  = x(K.expc.idx.x3)';

if ~isempty(x1)
    
    %scale vars to dual: (see p. 134 in Chares)
    x1t = x1;
    x1  = -x3;
    x2  = exp(1)*x2;
    x3  = -x1t;

    [FEXP,GEXP,HEXP,feasEXP] = BarrFuncExp(x1,x2,x3,K,want);

    if feasEXP < 0
        F{4} = -1;
        return
    end
    
else
    FEXP = [ ]; GEXP = [ ]; HEXP = [ ];
end

%% (DUAL) Logarithmic cone variables

x1  = x(K.logc.idx.x1)';
x2  = x(K.logc.idx.x2)';
x3  = x(K.logc.idx.x3)';

if ~isempty(x1)
    
    %scale vars to dual: HOW??
    % THIS IS NOT CORRECT RIGHT NOW!
    %x1t = x1;
    %x1  = -x3;
    %x2  = exp(1)*x2;
    %x3  = -x1t;

    %XXX: Changed akle

    x1t = x1;
    x1  = -x2;
    x2  = exp(1)*x3;
    x3  = -x1t;



    [FLOG,GLOG,HLOG,feasLOG] = BarrFuncLog(x1,x2,x3,K,want);

    if feasLOG < 0
        F{4} = -1;
        return
    end
    
else
    FLOG = [ ]; GLOG = [ ]; HLOG = [ ];
end


%% Gather all

% assign output.
% if we got this far without returning, point is feasible
if want(1) > 0
    F{1} = accumarray(K.allsubsF,[FPOS;FPOW;FEXP;FLOG]);
end
if want(2) > 0
    F{2} = accumarray(K.allsubsG,[GPOS;GPOW;GEXP;GLOG]);
end
if want(3) > 0
    F{3} = accumarray(K.allsubsH,[HPOS;HPOW;HEXP;HLOG],[ ],[ ],[ ],true);
end
F{4} = 0;


function F = BarrFuncP(x,K,want)

% want(1): F{1} = Function value
% want(2): F{2} = Gradient
% want(3): F{3} = Hessian
%          F{4} = Feasibility (0 for feas, -1 for infeas)

%% Positive cone variables:

xpos = x(K.posc.idx.x);

if ~isempty(xpos)
    
    [FPOS,GPOS,HPOS,feasPOS] = BarrFuncPos(xpos,K,want);

    if feasPOS < 0
        F{4} = -1;
        return
    end
    
else
    FPOS = []; GPOS = []; HPOS = [];
end


%% Power cone variables:
x1  = x(K.powc.idx.x1)';
x2  = x(K.powc.idx.x2)';
x3  = x(K.powc.idx.x3)';

if ~isempty(x1)

    if K.powc.fpownest
        [FPOW,GPOW,HPOW,feasPOW] = BarrFuncPowNest(x1,x2,x3,K,want);
    elseif K.powc.fpowglin
        [FPOW,GPOW,HPOW,feasPOW] = BarrFuncPowGlin(x1,x2,x3,K,want);
    elseif K.powc.fpowglin2
        [FPOW,GPOW,HPOW,feasPOW] = BarrFuncPowGlin2(x1,x2,x3,K,want);
    elseif K.powc.fpowglin3
        [FPOW,GPOW,HPOW,feasPOW] = BarrFuncPowGlin3(x1,x2,x3,K,want);    
    elseif K.powc.fpowhild
        [FPOW,GPOW,HPOW,feasPOW] = BarrFuncPowRH(x1,x2,x3,K,want);  
    else
        [FPOW,GPOW,HPOW,feasPOW] = BarrFuncPow(x1,x2,x3,K,want);
    end
    

    if feasPOW < 0
        F{4} = -1;
        return
    end

else
    FPOW = []; GPOW = []; HPOW = [];
end

%% Exponential cone variables:
x1  = x(K.expc.idx.x1)';
x2  = x(K.expc.idx.x2)';
x3  = x(K.expc.idx.x3)';

if ~isempty(x1)

    [FEXP,GEXP,HEXP,feasEXP] = BarrFuncExp(x1,x2,x3,K,want);

    if feasEXP < 0
        F{4} = -1;
        return
    end
    
else
    FEXP = []; GEXP = []; HEXP = [];
end

%% Logarithmic cone variables
x1  = x(K.logc.idx.x1)';
x2  = x(K.logc.idx.x2)';
x3  = x(K.logc.idx.x3)';

if ~isempty(x1)

    [FLOG,GLOG,HLOG,feasLOG] = BarrFuncLog(x1,x2,x3,K,want);

    if feasLOG < 0
        F{4} = -1;
        return
    end
    
else
    FLOG = []; GLOG = []; HLOG = [];
end


%% Gather all:

% assign output.
% if we got this far without returning, point is feasible
if want(1) > 0
    F{1} = accumarray(K.allsubsF,[FPOS;FPOW;FEXP;FLOG]);
end
if want(2) > 0
    F{2} = accumarray(K.allsubsG,[GPOS;GPOW;GEXP;GLOG]);
end
if want(3) > 0
    F{3} = accumarray(K.allsubsH,[HPOS;HPOW;HEXP;HLOG],[],[],[],true);
end
F{4} = 0;







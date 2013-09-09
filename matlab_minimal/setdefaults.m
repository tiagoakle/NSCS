function pars = setdefaults(pars)

pars.ndefaultsused = 0;

% reals:
pars = setonedef('beta'       ,0.2,0,inf,'r',pars);
pars = setonedef('rhoP'       ,1e-6,0,inf,'r',pars);
pars = setonedef('rhoD'       ,1e-6,0,inf,'r',pars);
pars = setonedef('rhoA'       ,1e-6,0,inf,'r',pars);
pars = setonedef('rhoG'       ,1e-6,0,inf,'r',pars);
pars = setonedef('rhoI'       ,1e-6,0,1,'r',pars);
pars = setonedef('rhoM'       ,1e-8,0,inf,'r',pars);
pars = setonedef('rhoC1'      ,1e-8,0,inf,'r',pars);
pars = setonedef('rhoC2'      ,1e-8,0,inf,'r',pars);
pars = setonedef('theta'      ,0.80,0,inf,'r',pars);
pars = setonedef('eta'        ,0.9995,0.5,1.0,'r',pars);
pars = setonedef('lsccent'    ,0.50,0,inf,'r',pars);
pars = setonedef('lscaff'     ,0.94,0,inf,'r',pars);
pars = setonedef('stpszmin'   ,1e-6,0,1,'r',pars);
pars = setonedef('affstep1'   ,0.995,0,inf,'r',pars);
pars = setonedef('affstep2'   ,1.1,0,inf,'r',pars);
pars = setonedef('corstep1'   ,1.0,0,inf,'r',pars);
pars = setonedef('RK2param'   ,0.7,-0.01,1.01,'r',pars);
pars = setonedef('hlim'       ,0.9,-0.01,1.01,'r',pars);
pars = setonedef('checkDpt'   ,1,-2,2,'i',pars);
pars = setonedef('dualcent'   ,0,-1,3,'i',pars);

% integers:
pars = setonedef('outermaxit' ,500,1,inf,'i',pars);
pars = setonedef('innermaxit' ,50,1,inf,'i',pars);
pars = setonedef('lsmaxit'    ,300,1,inf,'i',pars);
pars = setonedef('nmaxexpand' ,300,1,inf,'i',pars);
pars = setonedef('nmaxsect'   ,300,1,inf,'i',pars);
pars = setonedef('nnup'       ,3,-0.5,50,'i',pars);
pars = setonedef('qntype'     ,1,0,3,'i',pars);
pars = setonedef('lifttype'   ,2,0,3,'i',pars);
pars = setonedef('permuteM'   ,1,-1,2,'i',pars);
pars = setonedef('secord'     ,0,-1,3,'i',pars);
pars = setonedef('secord2'    ,0,-1,3,'i',pars);
pars = setonedef('bfgsstop'   ,1,-1,3,'i',pars);
pars = setonedef('addcent'    ,0,-1,3,'i',pars);


%testing:
pars = setonedef('secord3'    ,0,-1,3,'i',pars);
pars = setonedef('secord4'    ,0,-1,3,'i',pars);
pars = setonedef('cnbfgsstps' ,0,-1,999,'i',pars);
pars = setonedef('inclmu'     ,0,-1,3,'i',pars);
pars = setonedef('scalpo'     ,0,-1,3,'i',pars);
pars = setonedef('fpownest'   ,0,-1,3,'i',pars);
pars = setonedef('fpowglin'   ,0,-1,3,'i',pars);
pars = setonedef('fpowglin2'  ,0,-1,3,'i',pars);
pars = setonedef('fpowglin3'  ,0,-1,3,'i',pars);
pars = setonedef('fpowhild'   ,0,-1,3,'i',pars);
pars = setonedef('useamax '   ,0,-1,3,'i',pars);
pars = setonedef('beta2'      ,0.1,0,inf,'r',pars);
pars = setonedef('centmeastype',2,0,inf,'r',pars);
pars = setonedef('cholinc'    ,1,-1,2,'i',pars);

% ----

% output choices (integers):
pars = setonedef('echo'       ,1,-30,30,'i',pars);
pars = setonedef('trace'      ,1,0,3,'i',pars);

end

function pars = setonedef(name,dval,lb,ub,type,pars)
if isfield(pars,name)
    cv = eval(['pars.',name,';']);

    switch type
        case 'r'
            if isreal(cv)
                tc = true;
            else
                tc = false;
            end
        case 'i'
            if cv == round(cv)
                tc = true;
            else
                tc = false;
            end

        otherwise
            tc = true;
    end

    if ~tc || cv < lb || cv > ub
        error([name,...
            ' has illegal value. Specify nothing to use default value.']);
    else
        v = cv;
    end
else
    v = dval;
    pars.ndefaultsused = pars.ndefaultsused + 1;
    pars.defaultsused{pars.ndefaultsused} = name;
end
eval(['pars.',name,'=v;']);
end



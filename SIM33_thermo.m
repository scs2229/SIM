% SIM_33thermo is a zero-dimensional model which calculates ice and
% liquid hydrometeor numbers in six bins. Included microphysics are
% rime-splintering, aggregation, breakup upon collision, droplet
% coalescence, and droplet shattering. Moist thermodynamic variables are 
% tracked as in Korolev and Mazin 2003 Appendix B. Thermodynamic variables,
% i.e. updraft velocity and initial temperature, are additional inputs.
% Code Development: Sylvia Sullivan 2016 - 2017, GIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
% duration = simulation length [s]           
% step     = time step between microphysics and thermodynamics solns [s] 
% etabr    = breakup weighting [decimal form not percent]
% etaagg   = ice aggregation weighting [decimal form not percent]
% etaRSg   = rime splintering on small graupel weighting
% etaRSG   = rime splintering on large graupel weighting
% etaCOA   = droplet coalescence weighting
% etaSH    = droplet shattering weighting
% INPstop  = threshold INP number [m-3]
% uz       = updraft velocity [m s-1]
% Tstart   = initial temperature [K]
% tau      = 6-element vector of growth times from one class to proceeding, 
%            i.e. taui, taug, tauf, taud, taur, tauR (ice phase first) [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs: 
% Nice     = total ice hydrometeor number, i.e., Ni (1) + Ng (2) + NG (3)
%            over simulation duration, NOTE: UNITS ARE [L-1]
%            possibility to do the same thing with Ndrop
% cont     = contribution of different microphysics (nucleation, NU; breakup,
%            BR; rime-splintering, RS; and shattering, SH) to the ice crystal
%            number
% INPt     = cumulative number of nucleated INP, NOTE: UNITS ARE [L-1]
% tauNUC   = time point at which nucleation stops (equivalenty, time point
%            when all INP have nuclated), optional 

function [Nice,cont,INPt,tauNUC] = SIM_33thermo(etabr,etaagg,etaRSg,etaRSG,...
    etaCOA,etaSH,duration,step,INPstop,uz,Tstart,tau)
    
    % set physical constants
    g    = 9.81;    % gravitational acceleration [m s-2]
    R    = 8.314;   % gas constant [J mol-1 K-1]
    MWa  = 0.02897; % molecular mass of air [kg mol-1]
    MWw  = 0.01802; % molecular mass of water vapor [kg mol-1]
    Ra   = R/MWa;   % gas constant for air [J kg-1 K-1]
    Rv   = R/MWw;   % gas constant for water vapor [J kg-1 K-1]
    cp   = 1850;    % heat capacity (Rogers & Yau Ch2) [J kg-1 K-1]
    rhoi = 919;     % density of ice (assumed -15 C) [kg m-3]
    rhow = 998;     % density of supercooled water (assumed -20 C) [kg m-3]
    rhoa = 1.395;   % density of air (assumed - 15 C) [kg m-3]
    Nccn = 10^8;    % coefficient in Twomey CCN spec: continental case HUCM 
                    % (see COSMO src_twomom line 17060)
    kccn = 0.308;     % exponent in Twomey CCN spec: continental case HUCM 
    set1 = false;   % booleans to turn off hydrometeor number evolution
    set2 = false;   % once that particular hydrometeor is exhausted
    set3 = false;
    set4 = false;
    imm  = 0.9;     % what percentage of INP immediately nucleate ice
    coal2 = 0.005;   % what percentage of the large droplets coalesce
    
    % set depositional and condensational growth lags
    tauNUC = 72*60;    % nucleation time, overestimated initially
    taui = tau(1);
    taug = tau(2);
    tauf = tau(3);
    taud = tau(4);
    taur = tau(5);
    tauR = tau(6);
    lags = [taui, taui + taug, taui + taug + tauf, taud, taud + taur, ...
        taud + taur + tauR];  
    
    % set initial conditons for moist thermodynamics
    init(1)  = 68000;           % pressure [Pa]
    init(4)  = Tstart;          % temperature [K]
    init(5)  = 10^(-8);         % supersaturation wrt liquid water [-]
    init(6)  = 10^(-3);         % vapor mixing ratio [-]
    init(7)  = 10^(-6);         % small droplet radius [m]
    init(8)  = 5*10^(-6);       % ice crystal radius [m]
    init(9)  = 5*10^(-5);       % small graupel major axis [m] 
    init(10) = 20*10^(-5);      % large graupel major axis [m]
    init(11) = 12*10^(-6);      % medium droplet radius [m]
    init(12) = 25*10^(-6);      % large droplet radius [m]
    Ni       = 0;               % ice crystal number [m-3]
    Nd       = 0;               % small droplet number [m-3]
    
    % moist thermodynamic variables for the 1st soln of number tendencies
    P0  = init(1);
    T0  = init(4);
    sw  = init(5); 
    qv  = init(6);
    rw0 = init(7);
    ri0 = init(8); 
    ag0 = init(9); 
    aG0 = init(10);  
    rr0 = init(11); 
    rR0 = init(12); 
    
    % moist thermodynamic tendencies for the 1st soln of number tendencies
    dqidt   = 4*pi*rhoi/rhoa*Ni*ri0*sw*iceGrowth(T0,P0);
    dqwdt   = 4*pi*rhow/rhoa*Nd*rw0*sw*liqGrowth(T0,P0);
    
    % variables to determine the supersaturation maximum
    dswdt    = 0;
    tsmax    = 10000;   
    counter2 = 1;
    
    % set mixing ratio ICs, constrained by Ni0 ,Ng0, NG0, ri0, ag0, and aG0
    mice    = rhoi*Ni*4/3*pi*ri0^3;
    mliq    = rhow*Nd*4/3*pi*rw0^3;
    init(2) = mliq;     %/(1 + mice); 
    init(3) = mice;     %/(1 + mliq);
    rhodep  = depDens(T0, P0, sw);
    
    % initiate the contribution and INP vectors
    cont     = struct('NU',0,'RSg',0,'RSG',0,'BR',0,'SH',0);
    INPt     = struct('x',zeros(1,15000),'y',zeros(1,15000));
    Nice     = zeros(1,15000);
    Ndrop    = zeros(1,15000);
    counter3 = 1;
    
    % set the simulation length and output structure
    timepts = linspace(0, duration, ceil(duration/step)); 
    soln    = struct('x',zeros(1,15000),'y',zeros(12,15000)); 
    counter = 1;
    
    clear sol  % ensure blank history
    % no graupel or medium or large droplets initially
    % solve microphysics for hydrometeor number [m-3]
    sol = dde23(@microphysics,lags,[Ni; 0; 0; Nd; 0; 0; T0; sw],[0,timepts(2)]);
    INPt.y(2) = cont.NU(end)*timepts(2);
        
    % save hydrometeor numbers [m-3] to pass to the thermodynamics
    Ni = sol.y(1,end); 
    Ng = sol.y(2,end); 
    NG = sol.y(3,end);
    Nd = sol.y(4,end); 
    Nr = sol.y(5,end); 
    NR = sol.y(6,end);
    
    % solve thermodynamics [all SI]
    sol2    = ode15s(@thermodynamics,[0,timepts(2)],init);
    l       = length(sol2.x);
    init(:) = sol2.y(:,end); 
        
    soln.x(1, counter:counter + l - 1) = sol2.x; 
    soln.y(:, counter:counter + l - 1) = sol2.y; 
    counter = counter + l;
        
    % save thermodynamic variables to pass to microphysics 
    P0  = init(1);
    T0  = init(4);
    sw  = init(5); sol.y(8,end) = sw;
    qv  = init(6);
    rw0 = init(7);
    ri0 = init(8); 
    ag0 = init(9); 
    aG0 = init(10);   
    rr0 = init(11); 
    rR0 = init(12);
        
    % iterate through the remaining time spans
    for tt = 2:length(timepts) - 1
        t0 = timepts(tt); tf = timepts(tt + 1);
        
        if (INPt.y(tt) >= INPstop && counter3 == 1)
            tauNUC   = tt; 
            counter3 = 2;
        end
        
        rhodep  = depDens(T0, P0, sw);
        options = ddeset('RelTol',1e-3,'AbsTol',1e-2); 
        sol = dde23(@microphysics,lags,sol,[t0,tf],options);
        INPt.x(tt+1) = tf;
        INPt.y(tt+1) = INPt.y(tt) + cont.NU(end)*(tf-t0);
        
        % save hydrometeor numbers [m-3] to pass to the thermodynamics
        Ni = sol.y(1,end);
        Ng = sol.y(2,end); 
        NG = sol.y(3,end);
        Nd = sol.y(4,end); 
        Nr = sol.y(5,end); 
        NR = sol.y(6,end);
        
        % large droplet coalescence at non-freezing temperatures
        if (etaSH ~= 0 && T0 > 273)
            NR  = (1 - coal2)*NR;
            rR0 = rR0*(1/(1 - coal2))^(1/3);
            sol.y(6,end) = NR;
            init(12) = rR0;
        end
        
        % calculate the total ice hydrometeor number from the sol structure
        % units are [L-1]
        Nice(tt+1) = (Ni + Ng + NG)./10^3;
        Ndrop(tt+1) = NR./10^3;
  
        % halt simulation if any hydrometeor numbers become negative
        if (NG < 0 || NR < 0) 
            break 
        end
        
        % reset if small or medium droplet, ice crystal or small graupel 
        % number become negative
        if (Nd < 0)
            Ndnow = sol.y(4,:);
            sol.y(4, Ndnow < 0) = 0;
            clear Ndnow indx
            Nd = 0;
            set1 = true;
        end
        
        if (Nr < 0)
            Nrnow = sol.y(5,:);
            sol.y(5, Nrnow < 0) = 0;
            clear Nrnow indx
            Nr = 0;
            set2 = true;
        end
        
        if (Ni < 0)
            Ninow = sol.y(1,:);
            sol.y(1, Ninow < 0) = 0;
            clear Ninow indx
            Ni = 0;
            set3 = true;
        end
        
        if (Ng < 0)
            Ngnow = sol.y(2,:);
            sol.y(2, Ngnow < 0) = 0;
            clear Ngnow indx
            Ng = 0;
            set4 = true;
        end
        
        % reset mixing ratios, constrained by new hydrometeor numbers
        mice = rhoi*Ni*4/3*pi*ri0^3 + rhodep*(Ng*4/3*pi*ag0^(iG(T0)+2) + ...
            NG*4/3*pi*aG0^(iG(T0)+2));
        mliq = rhow*4/3*pi*(Nd*rw0^3 + Nr*rr0^3 + NR*rR0^3);
        init(3) = mice; %/(1 + mliq); 
        init(2) = mliq; %/(1 + mice);
        
        % solve thermodynamics
        sol2    = ode15s(@thermodynamics,[t0,tf],init);
        l       = length(sol2.x);
        init(:) = sol2.y(:,end); 
        
        soln.x(1, counter:counter + l - 1) = sol2.x; 
        soln.y(:, counter:counter + l - 1) = sol2.y; 
        counter = counter + l;
        
        % save thermodynamic variables to pass to microphysics 
        T0  = init(4);
        sw  = init(5); sol.y(8,end) = sw;
        qv  = init(6);
        rw0 = init(7);
        ri0 = init(8); 
        ag0 = init(9); 
        aG0 = init(10);   
        rr0 = init(11); 
        rR0 = init(12);
        
        % reset the time to obtain max. supersaturation if dswdt < 0
        if (dswdt < 0 && counter2 == 1)
            tsmax   = tf;
            counter2 = 2;
        end
        
        % halt simulation if parcel becomes water subsaturated
        if (sw < 10^(-12) || T0 < 237)
            break
        end 
    end
    
    % remove trailing zeros from INPt
    indx   = find(INPt.y(1,:) ~= 0, 1, 'last'); 
    INPt.y = INPt.y(:,1:indx)./10^3; % units are [L-1]
    indx   = find(INPt.x(1,:) ~= 0, 1, 'last');
    INPt.x = INPt.x(:,1:indx); 
    indx   = find(Nice ~= 0, 1, 'last');
    Nice   = Nice(:,1:indx);
    indx   = find(Ndrop ~= 0, 1, 'last');
    Ndrop  = Ndrop(:,1:indx);
    
    % mixed phase supersaturation as in Korolev and Mazin 2003
    function dy = thermodynamics(~, y)
 
        % hydrostatic pressure evolution
        dy(1) = -g*uz*y(1)/(Ra*y(4));

        % droplet growth equation
        if (Nd > 0)
            dy(7) = y(5)/y(7)*liqGrowth(y(4),y(1));
        else
            dy(7) = 0;
        end
        
        % medium droplet growth equation
        if (Nr > 0)
            dy(11) = y(5)/y(11)*liqGrowth(y(4),y(1));
        else
            dy(11) = 0;
        end

        % large droplet growth equation
        if (NR > 0)
            dy(12) = y(5)/y(12)*liqGrowth(y(4),y(1));
        else
            dy(12) = 0;
        end

        % crystal growth equation
        if (ssICE(y(4),y(5)) > 0)
            dy(8) = ssICE(y(4),y(5))/y(8)*iceGrowth(y(4),y(1));
        else
            dy(8) = 0;
        end

        % small graupel growth equation
        if (Ng > 0 && ssICE(y(4),y(5)) > 0)
            dy(9) = ssICE(y(4),y(5))/axis(y(9),y(4))*graupGrowth(y(4),y(1),y(5));
        else
            dy(9) = 0;
        end

        % large graupel growth equation
        if (NG > 0 && ssICE(y(4),y(5)) > 0)
            dy(10) = ssICE(y(4),y(5))/axis(y(10),y(4))*graupGrowth(y(4),y(1),y(5));
        else
            dy(10) = 0;
        end

        % liquid mixing ratio evolution
        dy(2) = 4*pi*rhow/rhoa*(Nd*y(7)^2*dy(7) + Nr*y(11)^2*dy(11) + ...
            NR*y(12)^2*dy(12));

        % crystal mixing ratio evolution
        dy(3) = 4*pi/rhoa*(rhoi*Ni*y(8)^2*dy(8) + rhodep*Ng*axis(y(9),y(4))^2*dy(9) + ...
            rhodep*NG*axis(y(10),y(4))^2*dy(10));

        % temperature evolution
        dy(4) = -g*uz/cp + heatVap(y(4))/((1+y(6))*cp)*dy(2) + ...
            heatSub(y(4))/((1+y(6))*cp)*dy(3);

        % liquid supersaturation evolution
        dy(5) = (1 + y(5))*(g*uz*heatVap(y(4))/(cp*Rv*y(4)^2) - g*uz/(Ra*y(4)) - ...
            (1/y(6) + heatVap(y(4))^2/(cp*Rv*y(4)^2))*dy(2) - ...
            (1/y(6) + heatSub(y(4))*heatVap(y(4))/(cp*Rv*y(4)^2))*dy(3));
        dswdt = dy(5);
        
        % water vapor mixing ratio evolution
        dy(6) = -dy(2) - dy(3);
        
        dy = dy';  
    end

    function dy = microphysics(t, y, Z)
        % fragments generated upon collision
        Ncoll = Neject(y(7));
        % splinters generated, Hallett & Mossop 1974
        NRS = Nsplin(rR0);
        % fragments generated upon droplet shattering
        NSH = Nshatter(2*rR0);
        
        if (iG(y(7)) < 1) % oblate          
            xig = (1 - ag0^(1-iG(y(7))))*rhodep/rhoi + ag0^(1-iG(y(7))); 
        else % prolate
            xig = (1 - ag0^(iG(y(7))-1))*rhodep/rhoi + ag0^(y(7)-1); 
        end
        
        if (iG(y(7)) < 1)
            xiG = (1 - aG0^(1-iG(y(7))))*rhodep/rhoi + aG0^(1-iG(y(7))); 
        else
            xiG = (1 - aG0^(iG(y(7))-1))*rhodep/rhoi + aG0^(iG(y(7))-1);
        end
        
        % calculate terminal velocities for the larger hydrometeor classes
        vtr = termVel(rr0, y(7), rhodep); 
        vtR = termVel(rR0, y(7), rhodep); 
        vtg = termVel(ag0, y(7), rhodep); 
        vtG = termVel(aG0, y(7), rhodep); 
        
        % sweep-out kernels for rime-splintering [m3 s-1]      
        HMeta = HMeffic(y(7)); % fraction of cloud in the RS temperature zone
        KRSg  = pi*(xig*ag0^2 + rR0^2)*(vtg - vtR)*HMeta;
        KRSG  = pi*(xiG*aG0^2 + rR0^2)*(vtG - vtR)*HMeta; 
        
        % sweep-out kernel for collisions [m3 s-1]
        Kcoll = pi*(xig*ag0^2 + xiG*aG0^2)*vtG;  % - vtg
        
        % sweep-out kernel for aggregation of crystals and small graupel
        Kagg  = pi*(xig*ag0^2 + ri0^2)*vtg;
        
        % sweep-out kernel for droplet coalescence
        Kcoal = pi*(rr0^2 + rw0^2)*vtr;
        
        % probability of a droplet freezing and shattering
        pSH   = pShatter(y(7));
        pFR   = pFreeze(y(7), 2*rR0); 
        
        % hydrometeor numbers at the various delay
        ylag1 = Z(:,1); 
        ylag2 = Z(:,2); 
        ylag3 = Z(:,3);
        ylag4 = Z(:,4); 
        ylag5 = Z(:,5); 
        ylag6 = Z(:,6);
        
        % ice crystal number evolution; generation by collisions and
        % rime-splintering, loss to aggregation
        if (set3 == false)
            dy(1) = nucl(y(7))*imm*heaviside(tauNUC-t) + ...
                etabr*Kcoll*Ncoll*y(2)*y(3) + ...
                etaRSg*KRSg*NRS*y(2)*y(6) + ...
                etaRSG*KRSG*NRS*y(3)*y(6) + ...
                etaSH*pSH*pFR*NSH*y(6) - ...
                nucl(ylag1(7))*imm*heaviside(t-taui)*heaviside(tauNUC-t+taui) - ...
                etabr*Kcoll*Ncoll*ylag1(2)*ylag1(3) - ...
                etaRSg*KRSg*NRS*ylag1(2)*ylag1(6) - ...
                etaRSG*KRSG*NRS*ylag1(3)*ylag1(6) - ...
                etaSH*pSH*pFR*NSH*ylag1(6) - ...
                etaagg*Kagg*y(1)*y(2);     
        else
            dy(1) = 0;
        end
        
        % small graupel number evolution
        if (set4 == false)
            dy(2) = nucl(ylag1(7))*imm*heaviside(t-taui)*heaviside(tauNUC-t+taui) + ... 
                etabr*Kcoll*Ncoll*ylag1(2)*ylag1(3) + ...
                etaRSg*KRSg*NRS*ylag1(2)*ylag1(6) + ...
                etaRSG*KRSG*NRS*ylag1(3)*ylag1(6) + ...
                etaSH*pSH*pFR*NSH*ylag1(6) - ...
                nucl(ylag2(7))*imm*heaviside(t-taui-taug)*heaviside(tauNUC-t+taui+taug) - ... 
                etabr*Kcoll*Ncoll*ylag2(2)*ylag2(3) - ...
                etaRSg*KRSg*NRS*ylag2(2)*ylag2(6) - ...    
                etaRSG*KRSG*NRS*ylag2(3)*ylag2(6) - ...
                etaSH*pSH*pFR*NSH*ylag2(6) - ...
                etaagg*Kagg*y(1)*y(2); 
        else
            dy(2) = 0;
        end

        dy(3) = nucl(ylag2(7))*imm*heaviside(t-taui-taug)*heaviside(tauNUC-t+taui+taug) + ... 
            etabr*Kcoll*Ncoll*ylag2(2)*ylag2(3) + ...
            etaRSg*KRSg*NRS*ylag2(2)*ylag2(6) + ...
            etaRSG*KRSG*NRS*ylag2(3)*ylag2(6) + ...
            etaSH*pSH*pFR*NSH*ylag2(6) + ... 
            etaagg*Kagg*y(1)*y(2) - ...
            nucl(ylag3(7))*imm*heaviside(t-taui-taug-tauf)*heaviside(tauNUC-t+taui+taug+tauf) - ... 
            etabr*Kcoll*Ncoll*ylag3(2)*ylag3(3) - ...
            etaRSg*KRSg*NRS*ylag3(2)*ylag3(6) - ...
            etaRSG*KRSG*NRS*ylag3(3)*ylag3(6) - ...
            etaSH*pSH*pFR*NSH*ylag3(6);                        
        
        if (set1 == false)
            dy(4) = activ(y(7),y(8),t) - etaCOA*Kcoal*y(4)*y(5) - ...
                activ(ylag4(7),ylag4(8),t - taud)*heaviside(t - taud);
        else
            dy(4) = 0;
        end
        
        if (set2 == false)
            dy(5) = activ(ylag4(7),ylag4(8),t - taud)*heaviside(t - taud) - ...
                etaCOA*Kcoal*y(4)*y(5) - ...
                activ(ylag5(7),ylag5(8),t - taud - taur)*heaviside(t - taud - taur);
        else
            dy(5) = 0;
        end

        dy(6) = activ(ylag5(7),ylag5(8),t - taud - taur)*heaviside(t - taud - taur) + ...
            etaCOA*Kcoal*y(4)*y(5) - ...
            0.1*activ(ylag6(7),ylag6(8),t - taud - taur - tauR)*heaviside(t - taud - taur - tauR) - ...
            etaRSg*KRSg*y(2)*y(6) - ...
            etaRSG*KRSG*y(3)*y(6) - ...
            etaSH*pSH*pFR*y(6);
        
        dqidt = 4*pi/rhoa*ssICE(y(7),y(8))*iceGrowth(y(7),P0)*(rhoi*y(1)*ri0 + ...
            rhodep*y(2)*ag0 + rhodep*y(3)*aG0);
        dqwdt = 4*pi*rhow/rhoa*y(8)*liqGrowth(y(7),P0)*(y(4)*rw0 + ...
            y(5)*rr0 + y(6)*rR0);
        
        % temperature evolution
        % duplicated from thermodynamics for activation and nucleation 
        dy(7) = -g*uz/cp + heatVap(y(7))/((1 + qv)*cp)*dqwdt + ...
            heatSub(y(7))/((1 + qv)*cp)*dqidt;

        % liquid supersaturation evolution
        % duplicated from thermodynamics for activation and nucleation 
        dy(8) = (1 + y(8))*(g*uz*heatVap(y(7))/(cp*Rv*y(7)^2) - g*uz/(Ra*y(7)) - ...
            (1/qv + heatVap(y(7))^2/(cp*Rv*y(7)^2))*dqwdt - ...
            (1/qv + heatSub(y(7))*heatVap(y(7))/(cp*Rv*y(7)^2))*dqidt);

        % calculate contributions from different processes
        cont.NU = [cont.NU nucl(T0)*imm*heaviside(tauNUC-t)]; 
        cont.RSg = [cont.RSg etaRSg*KRSg*NRS*y(2)*y(6)]; 
        cont.RSG = [cont.RSG etaRSG*KRSG*NRS*y(3)*y(6)];
        cont.BR = [cont.BR etabr*Kcoll*Ncoll*y(2)*y(3)]; 
        cont.SH = [cont.SH etaSH*NSH*pSH*pFR*y(6)];
        
        dy = dy';
    end
    
    % correlation for heat of vaporization, Murphy and Koop Eq 9
    % 236 K < T < 273 K, Lw [=] J kg-1
    function Lw = heatVap(T)
        A(1) = 56579; A(2) = -42.212; A(3) = 0.1149; A(4) = 281.6;
        Lw = A(1) + A(2)*T + exp(A(3)*(A(4) - T));
        Lw = Lw/MWw;
    end

    % correlation for heat of sublimation, Rogers and Yau Ch2 
    % input T [=] C, Li [=] J kg-1
    function Li = heatSub(T)
        A(1) = 2834.1; A(2) = 0.29; A(3) = 0.004; 
        Tc = T - 273;
        Li = A(1) - A(2)*Tc - A(3)*Tc^2; Li = Li*10^3;  
    end
    
    % correlation for thermal conductivity of air, Kannuluik and Carman '51
    % input T [=] C, k [=] W m-1 K-1
    function k = kAir(T)
        Tc = T - 273;
        A(1) = 4.184*100; A(2) = 5.75*10^(-5); A(3) = 0.00317; 
        A(4) = 0.0000021;
        
        k = A(2)*(1 + A(3)*Tc - A(4)*Tc^2); 
        k = A(1)*k;         
    end

    % diffusivity of water vapor in air, Seinfeld and Pandis 17.61
    % D [=] m2 s-1
    function D = Diff(T, P)
        A(1) = 0.211; A(2) = 9.86923*10^(-6); A(3) = 1.94;
        D = A(1)/(P*A(2))*(T/273)^A(3);
        D = D/100^2;
    end

    % sat vapor pressure over water, Murphy and Koop
    % ew [=] Pa, previously used (www.srh.noaa.gov/)
    function ew = satVapW(T)
        A(1) = 54.842763; A(2) = -6763.22; A(3) = -4.21; A(4) = 0.000367;
        A(5) = 0.0415; A(6) = 218.8; A(7) = 53.878; A(8) = -1331.22; 
        A(9) = -9.44523; A(10) = 0.014025;
        ew = A(1) + A(2)/T + A(3)*log(T) + A(4)*T + ...
            atan(A(5)*(T - A(6)))*(A(7) + A(8)/T + A(9)*log(T) + A(10)*T);
        ew = exp(ew);
    end

    % temperature derivative of the sat vapor pressure over water
    function deweT = satVapWDeriv(T)
        A(2) = -6763.22; A(3) = -4.21; A(4) = 0.000367;
        A(5) = 0.0415; A(6) = 218.8; A(7) = 53.878; A(8) = -1331.22; 
        A(9) = -9.44523; A(10) = 0.014025;
        
        deweT = -A(2)/T^2 + A(3)/T + A(4) + A(5)/(A(5)^2*(T - A(6))^2 + 1)*...
            (A(7) + A(8)/T + A(9)*log(T) + A(10)*T) + atan(A(5)*(T - A(6)))*...
            (-A(8)/T^2 + A(9)/T + A(10));
    end

    % sat vapor pressure over ice, Murphy and Koop
    % ei [=] Pa, previously used (www.srh.noaa.gov/)
    function ei = satVapI(T)
        A(1) = 9.550426; A(2) = -5723.265; A(3) = 3.53068; A(4) = -0.00728332;
        ei = A(1) + A(2)/T + A(3)*log(T) + A(4)*T; 
        ei = exp(ei);
    end
    
    % calculate the ratio of sat vapor pressure over water to sat vapor
    % pressure over ice
    function si = ssICE(T, sw)
        xi = satVapW(T)/satVapI(T);
        si = xi*(sw + 1) - 1;
    end
    
    % calculate the growth factor for the liquid phase, Gw [=] m2 s-1
    function Gw = liqGrowth(T, P)
        Gw = rhow*heatVap(T)^2/(kAir(T)*Rv*T^2) + rhow*Rv*T/...
            (satVapW(T)*Diff(T,P));
        Gw = Gw^(-1);
    end
    
    % calculate the growth factor for the ice phase, Gi [=] m2 s-1
    function Gi = iceGrowth(T, P)
        Gi = rhoi*heatSub(T)^2/(kAir(T)*Rv*T^2) + rhoi*Rv*T/...
            (satVapI(T)*Diff(T,P));
        Gi = Gi^(-1);
    end
    
    % calculate the growth factor for the graupel, Gg [=] m2 s-1
    % deposition density correlation inside [=] kg m-3
    function Gg = graupGrowth(T, P, sw)
        rhodep = depDens(T,P,sw);
        
        Gg = rhodep*heatSub(T)^2/(kAir(T)*Rv*T^2) + ...
            rhodep*Rv*T/(satVapI(T)*Diff(T,P));
        Gg = Gg^(-1);
    end

    % calculate the deposition density
    function rhodep = depDens(T, P, sw)
        rhoexc = R*T*P*(1 - 1/(sw*satVapW(T)/satVapI(T))); 
        rhoexc = rhoexc/1000; % input in g m-3 to formula below
        
        rhodep = 0.91*exp(-3*max([rhoexc-0.05 0])/iG(T));
        rhodep = rhodep*1000; % convert to kg m-3 from g cc-1
    end

    % calculate the inherent growth ratio as a function of temperature,
    % polynomial fit to values from Chen and Lamb 1994
    function gamma = iG(T)
        if (T >= 243 && T < 254.75)
            a5 = 2.028163171467212*10^(-4);
            a4 = -0.2513674925340751;
            a3 = 124.6034934430437;
            a2 = -3.087994836362414*10^4;
            a1 = 3.826031142079379*10^6;
            a0 = -1.895991380396219*10^8;
            gamma = a5*T^5 + a4*T^4 + a3*T^3 + a2*T^2 + a1*T + a0;
        else
            if (T >= 254.75 && T < 268)
                a5 = -1.444706102657700e-04;
                a4 = 1.880196590316559e-01;
                a3 = -9.787081351890365e+01;
                a2 = 2.547063781212500e+04;
                a1 = -3.314084675273538e+06;
                a0 = 1.724705653261639e+08; 
                gamma = a5*T^5 + a4*T^4 + a3*T^3 + a2*T^2 + a1*T + a0;
            else
                if (T >= 268 && T < 273)
                    a4 = 2.975090226539706e-02;
                    a3 = -32.26566559477525;
                    a2 = 1.312224376894170e+04;
                    a1 = -2.371858356867212e+06;
                    a0 = 1.607667896119342e+08;
                    gamma = a4*T^4 + a3*T^3 + a2*T^2 + a1*T + a0;
                else
                    gamma = 1;
                end
            end
        end
    end

    % calculate the crystal / graupel capacitance as a function of the
    % hydrometeor's aspect ratio, McDonald 1963, c [=] input a
    function c = cap(a, T)
        if (iG(T) > 1) % prolate
            ee = sqrt(1 - a^(2*iG(T))/a^2);
            c = a*ee/asin(ee);
        else
            if (iG(T) < 1) % oblate
                AA = sqrt(a^(2*iG(T)) - a^2);
                c = AA/log((a^(iG(T)) + AA)/a);
            else 
                c = 1;
            end
        end
    end
    
    % use this to calculate a 'radius proxy' for spheroidal hydrometeors
    function rproxy = axis(a, T)
        c = cap(a, T);
        temp = (1 + iG(T));
        rproxy = (temp + 1)*a^(temp);
        rproxy = rproxy/(3*c);        
    end

    % calculate the portion of cloud which falls in the HM regime as in
    % Ferrier 1994, T [=] K
    function eta = HMeffic(T)
        if (T <= 271 && T >= 269)
            HMT = 0.5;
        else
            if (T < 269 && T >= 267)
                HMT = 1;
            else
                if (T < 267 && T >= 265)
                    HMT = 0.5;
                else
                    HMT = 0.05;
                end
            end
        end

        eta = HMT;
    end

    % determine the number of ice fragments ejected upon collision
    % according to a four-parameter function (Takahashi et al. 1995)
    function frag = Neject(T)
        F = 280; Tmin = 252; gamma = 5;
        frag = F*(T - Tmin)^(1.2)*exp(-(T - Tmin)/gamma);
        frag(frag < 0) = 0; %10
    end

    % determine the number of ice fragments splintered off per number of
    % large droplets as a function of those rimable droplets' diameter
    % constant from Hallett and Mossop 1974, per mg-1 of rime
    function frag = Nsplin(D)
        frag = 3*10^8*rhow*pi/6*D^3;
    end

    % determine the number of ice fragments ejected upon shattering of a
    % frozen droplet as a function of droplet size (Lawson et al. 2014)
    % D [=] m
    function frag = Nshatter(D)
        D = D/10^(-6);
        frag = 2.5*10^(-11)*D^4;
    end
    
    % probability of a liquid droplet shattering during freezing (based
    % approximately on A. Kiselev experiments), T [=] K
    function p = pShatter(T)
        if (T < 266)
            mu  = 258; sigma = 5;
            max = 0.2; % normalize to a max probability of 20%
            p   = normpdf(T,mu,sigma).*max;
        else
            p = 0;
        end
    end

    % probability of a liquid droplet freezing, Paukert et al. 2017
    % T [=] K, droplet diameter D [=] m
    function pf = pFreeze(T, D)
       if (Nd ~= 0 && D > 100*10^(-6) && T < 273)
           f    = 1;
           inpT = (0.117*exp(-0.125*(T - 273.2))/f)*10^3;
                % above 10^3 converts L-1 to m-3
           pf   = (1 - imm)*inpT/Nd*50;
       else
           pf = 0;
       end
    end

    % calculate the terminal velocity of a hydrometer, based on Mitchell
    % and Heymsfield 2005, function of temperature and maximum dimension
    function vt = termVel(dim, T, rhod)

        % calculate the Best (Davies) number as in Lamb and Verlinde(9.29)
        mu = viscosity(T);
        X = 8*rhoa*rhod*g/(3*mu^2)*dim^3;
        C0 = 0.6; delta0 = 5.83; % as in Khvorostyanov and Curry 2002
        term = 4*sqrt(X)/(delta0^2*sqrt(C0));
        Re = delta0^2/4*(sqrt(1 + term) - 1)^2;
            % a0X^b0 term of MH05 is neglected for graupel and hail
        vt = Re/dim*mu/rhoa;
    end

    % calculate the viscosity of air as a function of temperature according
    % to the Sutherland model, [Pa-s]
    function mu = viscosity(T)
        mu0 = 18.27*10^(-6); % [Pa-s]
        C = 120; To = 291.15; % [K]
        mu = mu0*(To + C)/(T + C)*(T/To)^(3/2);
    end
    
   % droplet activation rate based on the Twomey cont. CCN spec [m-3 s-1]
    function newNd = activ(T, sw, tiempo)
        % below about -35 C homogeneous nucleation can occur
        % beyond the supersaturation maximum no activation occurs
        if (T > 250 && tiempo <= tsmax)  
            lapse    = -6/1000;
                % moist adiabatic lapse rate of 6 K/km -> [K m-1]
            dswdtEST = sw*satVapWDeriv(T);
            dswdtEST = dswdtEST + (satVapW(T)*lapse)^(-1)*g*rhoa;
            dswdtEST = dswdtEST*lapse*uz;        
            newNd    = Nccn*kccn*(sw*100)^(kccn-1)*dswdtEST;
        else
            newNd = 0;
        end
    end
    
    % ice nucleation rate based on the DeMott et al. 2010 spectrum
    function newNi = nucl(T0)
        f = 1;
        if (T0 > 237 && T0 < 273)
            aa = -14.625; 
                % -14.625 = 0.117*-0.125*1000 [m-3 K-1]
            lapse = -6/1000;
                % moist adiabatic lapse rate of 6 K/km -> [K m-1]
            newNi = aa*uz*lapse*exp(-0.125*(T0 - 273.2))/f;
        else 
            newNi = 0;
        end
    end
end

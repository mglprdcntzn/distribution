%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3 for distribution
% Miguel Parada Contzen
% Conce, enero 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc
rng(1); %random numbers seed
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grafgrabar   = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%some drawing constants
    ancho = 200;
    alto  = 200;
    Fontsize = 24;
    LineSize = 2;

    hrs  = 0:24;
    hrtext  = {'00:00','','','','','',...
               '06:00','','','','','',...
               '12:00','','','','','',...
               '18:00','','','','','','24:00'} ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general defitions
    j = sqrt(-1);
    V = 15e3;%4.16e3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nodes
    NodeData = dlmread('nodes.csv',';',1,0);
    N = size(NodeData,1) - 1;

    unos = ones(N,1);
    Id   = eye(N,N);
%Nominal load
    Ploadmax  = NodeData(2:end,4)*1e3;
%Load mix
    LoadMix = NodeData(2:end,5:7)/100;
%RES
    Sresinstal = [NodeData(2:end,9),NodeData(2:end,8)]*1e3; %Wind; PV
%Compensation
    Qcomp = NodeData(2:end,10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Branches
    BranchData = dlmread('branches.csv',';',1,0);
    NB = size(BranchData,1);
    
    DDD = zeros(NB,N+1); %incidenz matrix
    WWW = zeros(NB,NB); %branch admitances
    for bb = 1:NB,
        n = BranchData(bb,1)+1;
        m = BranchData(bb,2)+1;
        DDD(bb,n) =  1;
        DDD(bb,m) = -1;
        WWW(bb,bb) = 1/(BranchData(bb,3) + j*BranchData(bb,4));
    end
%addmitances matrices
    hatYYY = DDD'*WWW*DDD;
    Y00    = hatYYY(1,1);
    Y0     = diag(hatYYY(1,2:N+1));
    YYY    = hatYYY(2:end,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Profiles in p.u.
    Phr = dlmread('LoadProfiles.csv',';',1,1);
    Phr = [Phr;Phr(1,:)]; %00:00 = 24:00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm tunning
    tol        = 1e-6; %tolerance
    maxit      = 100; %max number iterations
    terminator = 0;%0: h_I; 1: h_P
        fprintf('   Max. Number of iterations: %i\n',maxit)
        fprintf('   Tolerance: %e\n',tol)

    algoritmo = 'CCSS';
        
    invMtx = YYY\eye(size(YYY));
    bbb    = invMtx*Y0*unos*V;   
        fprintf('   Solution algorithm: %s (%1.0f)\n',algoritmo,5);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exponential function (expm too slow?)
    my_expmj = @(XXX) diag(  cos(diag(XXX)) + j*sin(diag(XXX)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power balance functions
    h_P = @(VVV,Sres,Sz,Si,Sp)   Sres*unos -Sz*(abs(VVV)^2)*unos - Si*abs(VVV)*unos - Sp*unos - V*conj(Y0)*VVV*unos - VVV*conj(YYY)*conj(VVV)*unos;
    A_P = @(VVV,Sres,Sz,Si,Sp)  -Sz*conj(VVV)-0.5*Si*conj(VVV)/abs(VVV) - V*conj(Y0) - diag(conj(YYY)*conj(VVV)*unos);
    B_P = @(VVV,Sres,Sz,Si,Sp)  -Sz*VVV - 0.5*Si*VVV/abs(VVV)-VVV*conj(YYY);
    C_P = @(AAA,BBB,RRR,Psi)  AAA*my_expmj(Psi) + BBB*my_expmj(-Psi);
    D_P = @(AAA,BBB,RRR,Psi)  AAA*j*RRR*my_expmj(Psi) - BBB*j*RRR*my_expmj(-Psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current balance functions
    h_I = @(VVV,Sres,Sz,Si,Sp)  VVV\Sres*unos -Sz*conj(VVV)*unos - (Si*abs(VVV)/VVV)*unos- VVV\Sp*unos - conj(Y0)*unos*V - conj(YYY)*conj(VVV)*unos;
    A_I = @(VVV,Sres,Sz,Si,Sp)  -(VVV^2)\Sres +Si*(conj(VVV))^2/(abs(VVV)^3) + (VVV^2)\Sp;
    B_I = @(VVV,Sres,Sz,Si,Sp)  -Sz - Si/abs(VVV)-conj(YYY);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terminator functions
    if terminator,
        TheEnd = @(VVV,Sres,Sz,Si,Sp) norm(h_P(VVV,Sres,Sz,Si,Sp))<tol;
        fprintf('   Terminator: h_P\n')
    else
        TheEnd = @(VVV,Sres,Sz,Si,Sp) norm(h_I(VVV,Sres,Sz,Si,Sp))<tol;
        fprintf('   Terminator: h_I\n')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ZIP and fp values (Taken from Kundur:1994, Ch. 7)
    fpload = [ 0.85; %Com Summer
           0.90; %Com Winter
           0.90; %Res Summer
           0.99; %Res Winter
           0.85 %Ind
           ];
    aP     = [ 0.99; %Com Summer
           1.30; %Com Winter
           1.20; %Res Summer
           1.50; %Res Winter
           0.18 %Ind
           ];
    aQ     = [ 3.50; %Com Summer
           3.10; %Com Winter
           2.90; %Res Summer
           3.20; %Res Winter
           6.00 %Ind
           ];
    pz = 0.50; %arbitrary for all classes and nodes
    qz = 0.50; %arbitrary for all classes and nodes

    epsilon = 0.01; %for noise component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time definitions
    DT     = 10;%0.5;%1;% %sample period in minutes
    t0     = 0; %begining 00:00hrs
    tf     = 24*60; %end 24:00hrs in minutes
    t      = t0:DT:tf; % time in minutes
    
    % repeat initial instant for IC calculation
    reps   = 100;
    t      = [t0*ones(1,reps),t]; 
    kinic  = reps+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seasons of the year
    mes    = 10; %month of the year
    fprintf('   Month of the year: %i\n',mes)
    %Load parameters
    if mes>11 || mes <3, %if summer
        season = 'summer';
        
        fpcom = fpload(1);
        fpres = fpload(3);

        aPcom = aP(1);
        aPres = aP(3);
        aQcom = aQ(1);
        aQres = aQ(3);
        
        nuP_com = 1.0;
        nuP_res = 1.0;
    elseif mes <6, %if autumn
        season = 'autumn';
        
        fpcom = fpload(2);
        fpres = fpload(4);

        aPcom = aP(2);
        aPres = aP(4);
        aQcom = aQ(2);
        aQres = aQ(4);
        
        nuP_com = 1.2;
        nuP_res = 1.3;
    elseif mes <9, %if winter
        season = 'winter';
        
        fpcom = fpload(2);
        fpres = fpload(4);

        aPcom = aP(2);
        aPres = aP(4);
        aQcom = aQ(2);
        aQres = aQ(4);
        
        nuP_com = 1.4;
        nuP_res = 1.6;
    else %if spring
        season = 'spring';
        
        fpcom = fpload(1);
        fpres = fpload(3);

        aPcom = aP(1);
        aPres = aP(3);
        aQcom = aQ(1);
        aQres = aQ(3);
        
        nuP_com = 1.2;
        nuP_res = 1.3;
    end
        fpind = fpload(5);
        aPind = aP(5);
        aQind = aQ(5);
    fprintf('   Season: %s\n',season)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RES generation profiles
    %Photo-voltaic
    PerfilesSol      = csvread('PerfilSolCarrielSur.csv');
    Imax             = max(PerfilesSol(:,1));
    Perfiles         = zeros(25,12);
    for m=1:12,
        vec = PerfilesSol(((m-1)*24+1):(m*24),1);
        vec = [vec;vec(1)];
        Perfiles(:,m) = vec;
        
        if m == mes,
            PerfilMensualIrr = vec;
        end 
    end
    %Wind profile   
    PerfilesEol      = csvread('PerfilEolCarrielSur.csv');
    Wndmax           = max(PerfilesEol(:,2)).^3;
    Perfiles         = zeros(25,12);
    for m=1:12,
        vec = PerfilesEol(((m-1)*24+1):(m*24),2);
        vec = [vec;vec(1)];
        Perfiles(:,m) = vec;
        
        if m == mes,
            PerfilMensualWnd = vec;
        end 
    end
    PerfilMensualWnd = PerfilMensualWnd.^3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN LOOP for comparison
    fprintf('   Open Loop\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power factor control
    fpref    = 0.90;
    fpmin    = 0.80;
    Kfp      = 0*0.01;
    fprintf('   Proportional constant of fp ctrl: %1.2f\n',Kfp)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocating vectors
    Pz  = zeros(N,length(t));
    Qz  = zeros(N,length(t));
    Pi  = zeros(N,length(t));
    Qi  = zeros(N,length(t));
    Pp  = zeros(N,length(t));
    Qp  = zeros(N,length(t));

    Psolar   = zeros(N,length(t));
    fpsolar  = ones(N,length(t));
    Qsolar   = zeros(N,length(t));

    Peol     = zeros(N,length(t));
    fpeol    = ones(N,length(t));
    Qeol     = zeros(N,length(t));

    ve       = zeros(N,length(t));
    vant     = V*unos;

    P0       = zeros(1,length(t));
    Q0       = zeros(1,length(t));
    fp0      = zeros(1,length(t));

    SPfp     = zeros(1,length(t));
    SPfp(1)  = 1;%fpref;%
    m        = zeros(1,length(t));
    b        = zeros(1,length(t));
    sold     = 0*unos;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on time
for kk=1:length(t),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find time in hours
        hr  = floor(t(kk)/60);
        if hr>24,
            hr = hr - 24;
        end
        nxthr = min(hr+1,24);
        interhr = mod(t(kk),60)/60;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Define ZIP loads
        for ii = 1:N,%loop over nodes
            %comercial ZIP components
            ppucomi = (1-interhr)*Phr(hr+1,1) + interhr*Phr(nxthr+1,1);
            Pcomi = Ploadmax(ii)*LoadMix(ii,1)*ppucomi*(1+epsilon*randn)*nuP_com;
            Qcomi = Pcomi*sqrt((1/fpcom^2) - 1)*(1+epsilon*randn);

            if Pcomi>0,
                pIcomi = V*aPcom/Pcomi - 2*pz;
                qIcomi = V*aQcom/Qcomi - 2*qz;
            else
                pIcomi = 0;
                qIcomi = 0;
            end
            pPcomi = 1 - pz - pIcomi;
            qPcomi = 1 - qz - qIcomi;
            %residential ZIP components
            ppuresi = (1-interhr)*Phr(hr+1,2) + interhr*Phr(nxthr+1,2);
            Presi = Ploadmax(ii)*LoadMix(ii,2)*ppuresi*(1+epsilon*randn)*nuP_res;
            Qresi = Presi*sqrt((1/fpres^2) - 1)*(1+epsilon*randn);
            
            if Presi>0,
                pIresi = V*aPres/Presi - 2*pz;
                qIresi = V*aQres/Qresi - 2*qz;
            else
                pIresi = 0;
                qIresi = 0;
            end
            pPresi = 1 - pz - pIresi;
            qPresi = 1 - qz - qIresi;
            %industrial ZIP components
            ppuindi = (1-interhr)*Phr(hr+1,3) + interhr*Phr(nxthr+1,3);
            Pindi = Ploadmax(ii)*LoadMix(ii,3)*ppuindi*(1+epsilon*randn);
            Qindi = Pindi*sqrt((1/fpind^2) - 1)*(1+epsilon*randn);

            if Pindi>0,
                pIindi = V*aPind/Pindi - 2*pz;
                qIindi = V*aQind/Qindi - 2*qz;
            else
                pIindi = 0;
                qIindi = 0;
            end
            pPindi = 1 - pz - pIindi;
            qPindi = 1 - qz - qIindi;
            %Build load matrices
            Pz(ii,kk) = pz*(Pcomi + Presi + Pindi)/(V^2);
            Qz(ii,kk) = qz*(Qcomi + Qresi + Qindi)/(V^2);
            Pi(ii,kk) = (pIcomi*Pcomi + pIresi*Presi + pIindi*Pindi)/V;
            Qi(ii,kk) = (qIcomi*Qcomi + qIresi*Qresi + qIindi*Qindi)/V;
            Pp(ii,kk) = pPcomi*Pcomi + pPresi*Presi + pPindi*Pindi;
            Qp(ii,kk) = qPcomi*Qcomi + qPresi*Qresi + qPindi*Qindi;
        end
        %Complex ZIP loads
        Sz = diag(Pz(:,kk)+j*Qz(:,kk));
        Si = diag(Pi(:,kk)+j*Qi(:,kk));
        Sp = diag(Pp(:,kk)+j*Qp(:,kk));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Photo-voltaic Generation
        Irr = (1-interhr)*PerfilMensualIrr(hr+1) + interhr*PerfilMensualIrr(nxthr+1);
        Psolar(:,kk)  = (unos+epsilon*randn(N,1)).*Sresinstal(:,2)*Irr/Imax;
        fpsolar(:,kk) = SPfp(kk)*unos;
        Qsolar(:,kk)  = Psolar(:,kk).*sqrt((1./fpsolar(:,kk).^2) - 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Wind Power Generation
        Wnd = (1-interhr)*PerfilMensualWnd(hr+1) + interhr*PerfilMensualWnd(nxthr+1);
        Peol(:,kk)  = (unos+epsilon*randn(N,1)).*Sresinstal(:,1)*Wnd/Wndmax;
        fpeol(:,kk) = SPfp(kk)*unos;
        Qeol(:,kk)  = Peol(:,kk).*sqrt((1./fpeol(:,kk).^2) - 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RES complex matrix
        Sres = diag(j*Qcomp ... %compensation
                  + Psolar(:,kk) + j*Qsolar(:,kk) ... %PV
                  + Peol(:,kk) + j*Qeol(:,kk) ); %wind
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Initial condition for iterations
        VV   = diag(vant);
        RR   = abs(VV);
        PSI  = angle(VV);
        rp   = [RR*unos;PSI*unos];
        vvv  = VV*unos;
        uw   = [real(VV*unos);imag(VV*unos)];
            permiso = 1;
            cont    = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sucesive currents iteration loop
            while permiso,
                sss   = Sres*unos - Sz*(abs(vvv).^2) - Si*abs(vvv) - Sp*unos;
                iii   = conj(sss./vvv);
                vvv   = invMtx*iii - bbb ;
                VV    = diag(vvv);

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if TheEnd(VV,Sres,Sz,Si,Sp),
                    permiso = 0;
                end
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save calculated voltage
        ve(:,kk) = VV*unos;
        vant     = ve(:,kk);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find powers at substation root node 0
        s0      = conj(Y00)*V^2 + V*unos'*conj(Y0)*conj(vant);
        P0(kk)  = real(s0);
        Q0(kk)  = imag(s0);
        fp0(kk) = real(s0)/abs(s0);%*(1+epsilon*randn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Find slope "m" and free term "b" of the relation between pwr factors
        s0c       = conj(s0);
        alfa      = 1/abs(s0) - sqrt(s0c/(s0^3));
        H1        = alfa*conj(Y0) + conj(alfa)*Y0;
        H2        = -j*alfa*conj(Y0) + j*conj(alfa)*Y0;

        AA  = A_P(diag(ve(:,kk)),Sres,Sz,Si,Sp);
        BB  = B_P(diag(ve(:,kk)),Sres,Sz,Si,Sp);
        Mtx =  [real(AA+BB), -imag(AA-BB);
                imag(AA+BB),  real(AA-BB)];

        Pres = real(Sres);
        
        sss  = (1 + j*sqrt( (1/SPfp(kk)^2)- 1))*Pres*unos ...
               - abs(diag(ve(:,kk)))^2*Sz*unos ...
               - abs(diag(ve(:,kk)))*Si*unos - Sp*unos ;
        Dsss  = sss - sold;
        sold  = sss;

        m(kk) = (V/4)*unos'*[H1, H2]*inv(Mtx)*[0*unos;Pres*unos]./...
                ((SPfp(kk)^2)*sqrt(1 - (SPfp(kk)^2)));
        
        b(kk) = -(V/4)*unos'*[H1, H2]*inv(Mtx)*[real(Dsss);imag(Dsss)]/DT;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New power factor SP for RES
        e    = fpref-fp0(kk);
        temp = SPfp(kk) - DT*Kfp*e;
        SPfp(kk+1) = max(fpmin,min(temp,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Substation power factor
    figure(1)
        subplot(1,2,1)
            evalc('hold');
            plot([t(kinic),t(end)]/60, [fpref,fpref],':','LineWidth',LineSize*0.75)
            plot(t(kinic:end)/60,SPfp(kinic:end-1),'--','LineWidth',LineSize*0.75)
            plot(t(kinic:end)/60,fp0(kinic:end),'-','LineWidth',LineSize)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Open loop')
            legend('f_{p,ref}','f_{p,SP}','f_{p,0}','location','northeast')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLOSED LOOP
    fprintf('   Closed Loop\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power factor control
    Kfp      = 0.01;
    fprintf('   Proportional constant of fp ctrl: %1.2f\n',Kfp)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocating vectors
    Pz  = zeros(N,length(t));
    Qz  = zeros(N,length(t));
    Pi  = zeros(N,length(t));
    Qi  = zeros(N,length(t));
    Pp  = zeros(N,length(t));
    Qp  = zeros(N,length(t));

    Psolar   = zeros(N,length(t));
    fpsolar  = ones(N,length(t));
    Qsolar   = zeros(N,length(t));

    Peol     = zeros(N,length(t));
    fpeol    = ones(N,length(t));
    Qeol     = zeros(N,length(t));

    ve       = zeros(N,length(t));
    vant     = V*unos;

    P0       = zeros(1,length(t));
    Q0       = zeros(1,length(t));
    fp0      = zeros(1,length(t));

    SPfp     = zeros(1,length(t));
    SPfp(1)  = 1;%fpref;%
    m        = zeros(1,length(t));
    b        = zeros(1,length(t));
    sold     = 0*unos;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on time
for kk=1:length(t),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find time in hours
        hr  = floor(t(kk)/60);
        if hr>24,
            hr = hr - 24;
        end
        nxthr = min(hr+1,24);
        interhr = mod(t(kk),60)/60;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Define ZIP loads
        for ii = 1:N,%loop over nodes
            %comercial ZIP components
            ppucomi = (1-interhr)*Phr(hr+1,1) + interhr*Phr(nxthr+1,1);
            Pcomi = Ploadmax(ii)*LoadMix(ii,1)*ppucomi*(1+epsilon*randn)*nuP_com;
            Qcomi = Pcomi*sqrt((1/fpcom^2) - 1)*(1+epsilon*randn);

            if Pcomi>0,
                pIcomi = V*aPcom/Pcomi - 2*pz;
                qIcomi = V*aQcom/Qcomi - 2*qz;
            else
                pIcomi = 0;
                qIcomi = 0;
            end
            pPcomi = 1 - pz - pIcomi;
            qPcomi = 1 - qz - qIcomi;
            %residential ZIP components
            ppuresi = (1-interhr)*Phr(hr+1,2) + interhr*Phr(nxthr+1,2);
            Presi = Ploadmax(ii)*LoadMix(ii,2)*ppuresi*(1+epsilon*randn)*nuP_res;
            Qresi = Presi*sqrt((1/fpres^2) - 1)*(1+epsilon*randn);
            
            if Presi>0,
                pIresi = V*aPres/Presi - 2*pz;
                qIresi = V*aQres/Qresi - 2*qz;
            else
                pIresi = 0;
                qIresi = 0;
            end
            pPresi = 1 - pz - pIresi;
            qPresi = 1 - qz - qIresi;
            %industrial ZIP components
            ppuindi = (1-interhr)*Phr(hr+1,3) + interhr*Phr(nxthr+1,3);
            Pindi = Ploadmax(ii)*LoadMix(ii,3)*ppuindi*(1+epsilon*randn);
            Qindi = Pindi*sqrt((1/fpind^2) - 1)*(1+epsilon*randn);

            if Pindi>0,
                pIindi = V*aPind/Pindi - 2*pz;
                qIindi = V*aQind/Qindi - 2*qz;
            else
                pIindi = 0;
                qIindi = 0;
            end
            pPindi = 1 - pz - pIindi;
            qPindi = 1 - qz - qIindi;
            %Build load matrices
            Pz(ii,kk) = pz*(Pcomi + Presi + Pindi)/(V^2);
            Qz(ii,kk) = qz*(Qcomi + Qresi + Qindi)/(V^2);
            Pi(ii,kk) = (pIcomi*Pcomi + pIresi*Presi + pIindi*Pindi)/V;
            Qi(ii,kk) = (qIcomi*Qcomi + qIresi*Qresi + qIindi*Qindi)/V;
            Pp(ii,kk) = pPcomi*Pcomi + pPresi*Presi + pPindi*Pindi;
            Qp(ii,kk) = qPcomi*Qcomi + qPresi*Qresi + qPindi*Qindi;
        end
        %Complex ZIP loads
        Sz = diag(Pz(:,kk)+j*Qz(:,kk));
        Si = diag(Pi(:,kk)+j*Qi(:,kk));
        Sp = diag(Pp(:,kk)+j*Qp(:,kk));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Photo-voltaic Generation
        Irr = (1-interhr)*PerfilMensualIrr(hr+1) + interhr*PerfilMensualIrr(nxthr+1);
        Psolar(:,kk)  = (unos+epsilon*randn(N,1)).*Sresinstal(:,2)*Irr/Imax;
        fpsolar(:,kk) = SPfp(kk)*unos;
        Qsolar(:,kk)  = Psolar(:,kk).*sqrt((1./fpsolar(:,kk).^2) - 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Wind Power Generation
        Wnd = (1-interhr)*PerfilMensualWnd(hr+1) + interhr*PerfilMensualWnd(nxthr+1);
        Peol(:,kk)  = (unos+epsilon*randn(N,1)).*Sresinstal(:,1)*Wnd/Wndmax;
        fpeol(:,kk) = SPfp(kk)*unos;
        Qeol(:,kk)  = Peol(:,kk).*sqrt((1./fpeol(:,kk).^2) - 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RES complex matrix
        Sres = diag(j*Qcomp ... %compensation
                  + Psolar(:,kk) + j*Qsolar(:,kk) ... %PV
                  + Peol(:,kk) + j*Qeol(:,kk) ); %wind
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Initial condition for iterations
        VV   = diag(vant);
        RR   = abs(VV);
        PSI  = angle(VV);
        rp   = [RR*unos;PSI*unos];
        vvv  = VV*unos;
        uw   = [real(VV*unos);imag(VV*unos)];
            permiso = 1;
            cont    = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sucesive currents iteration loop
            while permiso,
                sss   = Sres*unos - Sz*(abs(vvv).^2) - Si*abs(vvv) - Sp*unos;
                iii   = conj(sss./vvv);
                vvv   = invMtx*iii - bbb ;
                VV    = diag(vvv);

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if TheEnd(VV,Sres,Sz,Si,Sp),
                    permiso = 0;
                end
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save calculated voltage
        ve(:,kk) = VV*unos;
        vant     = ve(:,kk);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find powers at substation root node 0
        s0      = conj(Y00)*V^2 + V*unos'*conj(Y0)*conj(vant);
        P0(kk)  = real(s0);
        Q0(kk)  = imag(s0);
        fp0(kk) = real(s0)/abs(s0);%*(1+epsilon*randn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Find slope "m" and free term "b" of the relation between pwr factors
        s0c       = conj(s0);
        alfa      = 1/abs(s0) - sqrt(s0c/(s0^3));
        H1        = alfa*conj(Y0) + conj(alfa)*Y0;
        H2        = -j*alfa*conj(Y0) + j*conj(alfa)*Y0;

        AA  = A_P(diag(ve(:,kk)),Sres,Sz,Si,Sp);
        BB  = B_P(diag(ve(:,kk)),Sres,Sz,Si,Sp);
        Mtx =  [real(AA+BB), -imag(AA-BB);
                imag(AA+BB),  real(AA-BB)];

        Pres = real(Sres);
        
        sss  = (1 + j*sqrt( (1/SPfp(kk)^2)- 1))*Pres*unos ...
               - abs(diag(ve(:,kk)))^2*Sz*unos ...
               - abs(diag(ve(:,kk)))*Si*unos - Sp*unos ;
        Dsss  = sss - sold;
        sold  = sss;

        m(kk) = (V/4)*unos'*[H1, H2]*inv(Mtx)*[0*unos;Pres*unos]./...
                ((SPfp(kk)^2)*sqrt(1 - (SPfp(kk)^2)));
        
        b(kk) = -(V/4)*unos'*[H1, H2]*inv(Mtx)*[real(Dsss);imag(Dsss)]/DT;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New power factor SP for RES
        e    = fpref-fp0(kk);
        temp = SPfp(kk) - DT*Kfp*e;
        SPfp(kk+1) = max(fpmin,min(temp,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Substation power factor
    figure(1)
        subplot(1,2,2)
            evalc('hold');
            plot([t(kinic),t(end)]/60, [fpref,fpref],':','LineWidth',LineSize*0.75)
            plot(t(kinic:end)/60,SPfp(kinic:end-1),'--','LineWidth',LineSize*0.75)
            plot(t(kinic:end)/60,fp0(kinic:end),'-','LineWidth',LineSize)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Closed loop')
            legend('f_{p,ref}','f_{p,SP}','f_{p,0}','location','northwest')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            ylim([fpmin,1]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % power factor slopes
    figure(2)
        subplot(1,2,1)
            plot(t(kinic:end)/60,m(kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            legend('m(\cdot)','Location','SouthEast')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,2,2)
            plot(t(kinic:end)/60,1+DT*Kfp*m(kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            legend('1+\Delta t  K_f m(\cdot)','Location','SouthEast')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save figures
if grafgrabar,
    sufix = ['_',season,'_mes_',num2str(mes),'.eps'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    figure(1) 
            set(gcf,'PaperUnits','centimeters',...
                    'PaperSize',[2*ancho 1*alto],...
                    'PaperPosition',[0 0 2*ancho 1*alto]); %[0 0 ancho alto]
            print('-depsc','-r200',['fp',sufix]) % FunciOn para guardar .eps 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
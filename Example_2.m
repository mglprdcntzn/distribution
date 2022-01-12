%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2 for distribution
% Miguel Parada Contzen
% Conce, enero 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc
rng(1); %random numbers seed
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ancho = 200;
    alto  = 200;
    Fontsize = 12;
    
    grafPerfiles = false;
    grafSimul    = false;
    grafgrabar   = false;
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
        
    algo = 1; %1: NR-ClAsico; 2: NR-UW; 3: NR-IUW; 4: NR-IUWsimp; 5: CCSS
    if algo==1,
        algoritmo = 'NR-classic';        
    elseif algo==2,
        algoritmo = 'NR-UW';       
    elseif algo==3,
        algoritmo = 'NR-IUW';
    elseif algo==4,
        algoritmo = 'NR-IUWsimp';
        
        Mtx = [-real(YYY), imag(YYY);
                imag(YYY), real(YYY)];
        invMtx = Mtx\eye(size(Mtx));
    elseif algo==5,
        algoritmo = 'CCSS';
        
        invMtx = YYY\eye(size(YYY));
        bbb    = invMtx*Y0*unos*V;
    end
        fprintf('   Solution algorithm: %s (%1.0f)\n',algoritmo,algo);
    
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
%Load profiles graph
    hrs  = 0:24;
    hrtext  = {'00:00','','','','','',...
               '06:00','','','','','',...
               '12:00','','','','','',...
               '18:00','','','','','','24:00'} ;
    if grafPerfiles,
        figure(1)
            subplot(1,3,1)
            evalc('hold');
            plot(hrs,Phr(:,1),'-','LineWidth',1.5)
            plot(hrs,Phr(:,2),'--','LineWidth',1.5)
            plot(hrs,Phr(:,3),'-.','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            axis([hrs(1),hrs(end),0.3,1.1])
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
            legend('Com.','Res.','Ind.',...
                   'location','south','orientation','horizontal') 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time definitions
    DT     = 0.5;%1;%10;% %sample period in minutes
    t0     = 0; %begining 00:00hrs
    tf     = 24*60; %end 24:00hrs in minutes
    t      = t0:DT:tf; % time in minutes
    
    % repeat initial instant for IC calculation
    reps   = 100;
    t      = [t0*ones(1,reps),t]; 
    kinic  = reps+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seasons of the year
    mes    = 7; %month of the year
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
    %Plot PV profiles
    if grafPerfiles,
        figure(1)
            subplot(1,3,2)
            plot(hrs,Perfiles,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('Mean Global radiation $[W/m^2]$','Interpreter','LaTeX','FontName','Times New Roman');
            legend('Jan','Feb','Mar',...
                   'Apr','May','Jun',...
                   'Jul','Agu','Sep',...
                   'Oct','Nov','Dec',...
                   'location','east',...
                   'orientation','vertical',...
                   'location','NorthWest')  
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
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
    %Plot wind profile
    if grafPerfiles,
        figure(1)
            subplot(1,3,3)
            hplot = plot(hrs,Perfiles,'-','LineWidth',1.5);
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('Mean wind speed $[m/s]$ at $10[m]$','Interpreter','LaTeX','FontName','Times New Roman');
            legend('Jan','Feb','Mar',...
                   'Apr','May','Jun',...
                   'Jul','Agu','Sep',...
                   'Oct','Nov','Dec',...
                   'location','east',...
                   'orientation','vertical',...
                   'location','NorthWest') 
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
            
            set(gcf,'PaperUnits','centimeters',...
                    'PaperSize',[3*ancho alto],...
                    'PaperPosition',[0 0 3*ancho alto]); %[0 0 ancho alto]
            if grafgrabar,
                print('-depsc','-r200','Perfiles.eps') % FunciOn para guardar .eps 
            end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power factor control
    fpref    = 0.90;
    fpmin    = 0.70;
    Kfp      = 0*0.1;
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

    IterationTime = zeros(1,length(t));
    Iterations    = zeros(1,length(t));
    
    SPfp     = zeros(1,length(t));
    SPfp(1)  = 1;%fpref;%
    m        = zeros(1,length(t));
    b        = zeros(1,length(t));
    SPfpant  = SPfp(1);
    fp0ant   = 1;
    sold     = 0*unos;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on time
beginsim = cputime;
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
        %select algorithm for solution
        if algo==1,
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NR-Classic iteration loop
            tic;
            while permiso,
                h   = h_P(VV,Sres,Sz,Si,Sp);

                AA  = A_P(VV,Sres,Sz,Si,Sp);
                BB  = B_P(VV,Sres,Sz,Si,Sp);
                CC  = C_P(AA,BB,RR,PSI);
                DD  = D_P(AA,BB,RR,PSI);

                Mtx =  [real(CC), real(DD);
                        imag(CC), imag(DD)];
                vtr = -[real(h);
                        imag(h)];

                Delta = Mtx\vtr;
                rp    = rp + Delta;
                RR    = diag(rp(1:N));
                PSI   = diag(rp(N+1:end));
                VV    = RR*my_expmj(PSI);

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if TheEnd(VV,Sres,Sz,Si,Sp),
                    permiso = 0;
                end
            end
        elseif algo==2,
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NR-UW iteration loop
            tic;
            while permiso,
                h   = h_P(VV,Sres,Sz,Si,Sp);

                AA  = A_P(VV,Sres,Sz,Si,Sp);
                BB  = B_P(VV,Sres,Sz,Si,Sp);

                Mtx =  [real(AA+BB), -imag(AA-BB);
                        imag(AA+BB),  real(AA-BB)];
                vtr = -[real(h);
                        imag(h)];

                Delta = Mtx\vtr;
                uw    = uw + Delta;
                VV    = diag(uw(1:N)) + j*diag(uw(N+1:end));

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if TheEnd(VV,Sres,Sz,Si,Sp),
                    permiso = 0;
                end
            end
        elseif algo==3,
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NR-IUW iteration loop
            tic;
            while permiso,
                h   = h_I(VV,Sres,Sz,Si,Sp);

                AA  = A_I(VV,Sres,Sz,Si,Sp);
                BB  = B_I(VV,Sres,Sz,Si,Sp);

                Mtx =  [real(AA+BB), -imag(AA-BB);
                        imag(AA+BB),  real(AA-BB)];
                vtr = -[real(h);
                        imag(h)];

                Delta = Mtx\vtr;
                uw    = uw + Delta;
                VV    = diag(uw(1:N)) + j*diag(uw(N+1:end));

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if TheEnd(VV,Sres,Sz,Si,Sp),
                    permiso = 0;
                end
            end
        elseif algo==4,
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NR-IUWsimpl iteration loop
            tic;
            while permiso,
                h   = h_I(VV,Sres,Sz,Si,Sp);

                vtr = -[real(h);
                        imag(h)];

                Delta = invMtx*vtr;
                uw    = uw + Delta;
                VV    = diag(uw(1:N)) + j*diag(uw(N+1:end));

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if TheEnd(VV,Sres,Sz,Si,Sp),
                    permiso = 0;
                end
            end
        elseif algo==5,
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sucesive currents iteration loop
            tic;
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
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time duration of algorithms
        elapsedTime = toc;
        IterationTime(kk) = elapsedTime;
        Iterations(kk) = cont;
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
endsim = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execution times
TotalTime = endsim - beginsim;
    fprintf('   Total simulation time: %1.2f [ms]\n',TotalTime*1000);
    fprintf('   Total iteration time: %1.2f [ms]\n',sum(IterationTime(kinic:end))*1000);
    % fprintf('   Time share ite/sim: %1.2f %%\n',100*sum(IterationTime)/TotalTime);
    fprintf('   Number of iterative processes: %i\n',length(kinic:length(t)));
    fprintf('   Iterations mean time: %1.2f [ms]\n',mean(IterationTime(kinic:end))*1000);
    fprintf('   Mean iterations: %1.2f\n',mean(Iterations(kinic:end)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post execution process
    %ZIP Loads
        Pload = Pz.*abs(ve).^2 + Pi.*abs(ve) + Pp;
        Qload = Qz.*abs(ve).^2 + Qi.*abs(ve) + Qp;
        fpload = Pload./sqrt(Qload.^2 + Pload.^2);
    %Power balance
        P  = Peol + Psolar - Pload;
        Q  = Qcomp*ones(1,length(t)) + Qeol + Qsolar - Qload;
        fp = P./sqrt(P.^2 + Q.^2);
        
        S     = P + j*Q;
        hP    = S - ve.*(V*conj(Y0)*unos*ones(1,length(t))) - ve.*conj(YYY*ve);
        hPr   = real(hP);
        hPi   = imag(hP);
        abshP = abs(hP);
        hP2   = hP.*conj(hP);
        
        hI    = (1./ve).*hP;
        hIr   = real(hI);
        hIi   = imag(hI);
        abshI = abs(hI);
        hI2   = hI.*conj(hI);
        
        fprintf('   Mean power balance error norm: %1.2e [W]\n',mean(sqrt(sum(hP2))));
        fprintf('   Mean current balance error norm: %1.2e [A]\n',mean(sqrt(sum(hI2))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Print table
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')   
    fprintf('  & %1.2f & %1.2f & %1.2e\n',mean(IterationTime(kinic:end))*1000,mean(Iterations(kinic:end)),mean(sqrt(sum(hI2))) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphs
if grafSimul,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sufix = ['_',season,'_mes_',num2str(mes),'_K_',num2str(Kfp),'_Alg_',algoritmo,'.eps'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ZIP components
    figure
        subplot(3,2,1)
            plot(t(kinic:end)/60,1000*Pz(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_Z [mS]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(3,2,2)
            plot(t(kinic:end)/60,1000*Qz(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Q_Z [mS]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(3,2,3)
            plot(t(kinic:end)/60,Pi(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_I [A]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(3,2,4)
            plot(t(kinic:end)/60,Qi(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Q_I [A]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(3,2,5)
            plot(t(kinic:end)/60,Pp(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_P [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(3,2,6)
            plot(t(kinic:end)/60,Qp(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Q_P [MVAr]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load powers
    figure
        subplot(1,2,1)
            plot(t(kinic:end)/60,Pload(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_{load} [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,2,2)
            plot(t(kinic:end)/60,Qload(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Q_{load} [MVAr]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
%         subplot(1,3,3)
%             plot(t(kinic:end)/60,fpload(:,kinic:end),'-','LineWidth',1.5)
%             box ON
%             set(gca,'XTick',hrs)
%             set(gca,'XTickLabel',hrtext)
%             title('fp_{load}')
%             grid ON
%                 hAx=gca;  % avoid repetitive function calls
%                 set(hAx,'xminorgrid','off','yminorgrid','off')
%             xlim([hrs(1),hrs(end)]);
%             set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PV generation
    figure
        subplot(1,3,1)
            plot(t(kinic:end)/60,Psolar(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_{PV} [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,2)
            plot(t(kinic:end)/60,Qsolar(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Q_{PV} [MVAr]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,3)
            plot(t(kinic:end)/60,fpsolar(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('f_{p,PV}')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Wind generation
    figure
        subplot(1,3,1)
            plot(t(kinic:end)/60,Peol(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_{wnd} [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,2)
            plot(t(kinic:end)/60,Qeol(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('Q_{wnd} [MVAr]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,3)
            plot(t(kinic:end)/60,fpeol(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('f_{p,wnd}')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RES alltogether
     figure
        subplot(1,2,1)
            plot(t(kinic:end)/60,Psolar(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_{PV} [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,2,2)
            plot(t(kinic:end)/60,Peol(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('P_{wnd} [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perturbations alltogether (Load-RES)
    figure
        subplot(1,3,1)
            plot(t(kinic:end)/60,P(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('real\{S\} [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,2)
            plot(t(kinic:end)/60,Q(:,kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('imag\{S\} [MVAr]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,3)
            plot(t(kinic:end)/60,fp(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('cos(\angle S)')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Voltages
    figure
        subplot(1,2,1)
            plot(t(kinic:end)/60,abs(ve(:,kinic:end))/V,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('|v| [pu]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,2,2)
            plot(t(kinic:end)/60,angle(ve(:,kinic:end))*180/pi,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('\angle v [°]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            set(gcf,'PaperUnits','centimeters',...
                    'PaperSize',[2*ancho 1*alto],...
                    'PaperPosition',[0 0 2*ancho 1*alto]); %[0 0 ancho alto]
            if grafgrabar,
                print('-depsc','-r200',['v',sufix]) % FunciOn para guardar .eps 
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Substation power
    figure
        subplot(1,2,1)
            plot(t(kinic:end)/60,P0(kinic:end)/1000000,'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('P_0 [MW]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
%         subplot(1,3,2)
%             plot(t(kinic:end)/60,Q0(kinic:end)/1000000,'-','LineWidth',1.5)
%             box ON
%             set(gca,'XTick',hrs)
%             set(gca,'XTickLabel',hrtext)
%             title('Q_0 [MVAr]')
%             grid ON
%                 hAx=gca;  % avoid repetitive function calls
%                 set(hAx,'xminorgrid','off','yminorgrid','off')
%             xlim([hrs(1),hrs(end)]);
%             set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,2,2)
%             evalc('hold');
            psi   = atan(Q0./P0);
            fpgen = cos(psi(kinic:end));%Pgen./sqrt(Pgen.^2 + Qgen.^2);
            tt    = t(kinic:end);
            plot(tt/60, fpgen,'-','LineWidth',1.0)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('f_{p,0}')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(gcf,'PaperUnits','centimeters',...
                    'PaperSize',[2*ancho 1*alto],...
                    'PaperPosition',[0 0 2*ancho 1*alto]); %[0 0 ancho alto]
            if grafgrabar,
                print('-depsc','-r200',['S0',sufix]) % FunciOn para guardar .eps 
            end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Power balance
    figure
        subplot(1,3,1)
            plot(t(kinic:end)/60,hPr(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('real\{h_{P}\} [W]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,2)
            plot(t(kinic:end)/60,hPi(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('imag\{h_{P}\} [VAr]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,3)
            plot(t(kinic:end)/60,abshP(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('|h_{P}| [VA]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Current balance
    figure
        subplot(1,3,1)
            plot(t(kinic:end)/60,hIr(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('real\{h_{I}\} [A]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,2)
            plot(t(kinic:end)/60,hIi(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('imag\{h_{I}\} [Ar]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,3)
            plot(t(kinic:end)/60,abshI(:,kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('|h_{I}| [|A|]')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Power factors relationship
    figure
        subplot(1,3,1)
            evalc('hold');
            plot([t(kinic),t(end)]/60, [fpref,fpref],':','LineWidth',1.0)
            plot(t(kinic:end)/60,SPfp(kinic:end-1),'--','LineWidth',1.0)
            plot(t(kinic:end)/60,fp0(kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('f_p')
            legend('f_{ref}','f_{SP}','f_{p,0}','location','northeast')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,2)
            evalc('hold');
            plot(t(kinic:end)/60,m(kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
%             title('m(\cdot)')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
        subplot(1,3,3)
            evalc('hold');
            plot(t(kinic:end)/60,1+DT*Kfp*m(kinic:end),'-','LineWidth',1.5)
            box ON
            set(gca,'XTick',hrs)
            set(gca,'XTickLabel',hrtext)
            title('1+\Delta t  K_f m(\cdot)')
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            xlim([hrs(1),hrs(end)]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
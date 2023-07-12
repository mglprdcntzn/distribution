%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for distribution with DG and H2
% Miguel Parada Contzen
% Conce, agosto 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc
rng(1); %random numbers seed
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grafgrabar   = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%some drawing constants
    ancho = 200;
    alto  = 200;
    Fontsize = 18;

    hrs  = 0:24;
    hrtext  = {'00:00','','','','','',...
               '06:00','','','','','',...
               '12:00','','','','','',...
               '18:00','','','','','','24:00'} ;
    hrtxtsh = {'00:00','','','','','',...
               ''     ,'','','','','',...
               '12:00','','','','','',...
               ''     ,'','','','','','24:00'} ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general defitions
    j = sqrt(-1);
    V = 15e3;%4.16e3;%
    
    fpmin    = 0.70;%min pwr fctr for DG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nodes
    NodeData = dlmread('nodesH2.csv',';',1,0);
    N = size(NodeData,1) - 1;

    unos = ones(N,1);
    Id   = eye(N,N);
%Nominal load
    Ploadmax  = NodeData(2:end,4)*1e3;
%Load mix
    LoadMix = NodeData(2:end,5:7)/100;
%RES
    SPVinstal = NodeData(2:end,8)*1e3; %PV
    PPVmax    = unos'*SPVinstal;
%Compensation
    Qcomp = NodeData(2:end,9)*1e3;
%H2
    PH2nom = NodeData(2:end,10)*1e3; %H2 storage
    PFCnom  = NodeData(2:end,11)*1e3; %H2 fuel cell
    
    PH2Tot  = sum(PH2nom);
    PFCTot  = sum(PFCnom);
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
    Ylines = WWW*DDD;
    hatYYY = DDD'*(WWW*DDD);
    Y00    = hatYYY(1,1);
    Y0     = diag(hatYYY(1,2:N+1));
    YYY    = hatYYY(2:end,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Profiles in p.u.
    Phr = dlmread('LoadProfiles.csv',';',1,1);
    Phr = [Phr;Phr(1,:)]; %00:00 = 24:00
    
    desfase = 30*(2*rand(N,1)-ones(N,1)); %15 mins forward or bck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power Flow Algorithm tunning
    tol        = 1e-6; %tolerance
    maxit      = 100; %max number iterations
        fprintf('   Max. Number of iterations: %i\n',maxit)
        fprintf('   Tolerance: %e\n',tol)

    invMtx = YYY\eye(size(YYY));
    bbb    = invMtx*Y0*unos*V;   
        fprintf('   Solution algorithm: %s (%1.0f)\n','CCSS',5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current balance functions
    h_I = @(VVV,Sres,Sp)  VVV\Sres*unos - VVV\Sp*unos - conj(Y0)*unos*V - conj(YYY)*conj(VVV)*unos;
    A_I = @(VVV,Sres,Sp)  -(VVV^2)\Sres + (VVV^2)\Sp;
    B_I = @(VVV,Sres,Sp)  -conj(YYY);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load fp values (Taken from Kundur:1994, Ch. 7)
    fpload = [ 0.85; %Com Summer
           0.90; %Com Winter
           0.90; %Res Summer
           0.99; %Res Winter
           0.85 %Ind
           ];

    epsilon = 0.01; %for noise component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time definitions
    h      = 5;%10;%0.5;% %sample period in minutes
    t0     = -12*60; %begining 00:00hrs
    tf     = 36*60; %end 24:00hrs in minutes
    t      = t0:h:tf; % time in minutes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seasons of the year
    mes    = 1; %month of the year / January=1 is summer
    fprintf('   Month of the year: %i\n',mes)
    %Load parameters
    if mes>11 || mes <3, %if summer
        season = 'summer';
        fpcom = fpload(1);
        fpres = fpload(3);
        
        nuP_com = 1.0;
        nuP_res = 1.0;
    elseif mes <6, %if autumn
        season = 'autumn';
        fpcom = fpload(2);
        fpres = fpload(4);
        
        nuP_com = 1.2;
        nuP_res = 1.3;
    elseif mes <9, %if winter
        season = 'winter';
        fpcom = fpload(2);
        fpres = fpload(4);
        
        nuP_com = 1.6;
        nuP_res = 1.7;
    else %if spring
        season = 'spring';
        fpcom = fpload(1);
        fpres = fpload(3);
        
        nuP_com = 1.2;
        nuP_res = 1.3;
    end
        fpind = fpload(5);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open Loop
    fprintf('   Open Loop\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocating vectors
    Pload  = zeros(N,length(t));
    Qload  = zeros(N,length(t));

    Psolar   = zeros(N,length(t));
    fpsolar  = ones(N,length(t));
    Qsolar   = zeros(N,length(t));

    ve       = zeros(N,length(t));
    vant     = V*unos;

    P0       = zeros(1,length(t));
    Q0       = zeros(1,length(t));
    fp0      = zeros(1,length(t));
    
    SPfp     = zeros(1,length(t));
    SPfp(1)  = 1;%fpref0;%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on time
for kk=1:length(t),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find time in hours
        hr  = floor((t(kk)+desfase)/60);
        hr  = mod(hr,24);
        if hr>24,
            hr = hr - 24;
        end
        nxthr = min(hr+1,24);
        interhr = mod((t(kk)+desfase),60)/60;
    %Define Power loads per class
        for ii = 1:N,%loop over nodes
            %comercial load components
            ppucomi = (1-interhr(ii))*Phr (hr(ii)+1,1) + interhr(ii)*Phr(nxthr(ii)+1,1);
            Pcomi = Ploadmax(ii)*LoadMix(ii,1)*ppucomi*(1+epsilon*randn)*nuP_com;
            Qcomi = Pcomi*sqrt((1/fpcom^2) - 1)*(1+epsilon*randn);
            %residential load components
            ppuresi = (1-interhr(ii))*Phr(hr(ii)+1,2) + interhr(ii)*Phr(nxthr(ii)+1,2);
            Presi = Ploadmax(ii)*LoadMix(ii,2)*ppuresi*(1+epsilon*randn)*nuP_res;
            Qresi = Presi*sqrt((1/fpres^2) - 1)*(1+epsilon*randn);
            %industrial load components
            ppuindi = (1-interhr(ii))*Phr(hr(ii)+1,3) + interhr(ii)*Phr(nxthr(ii)+1,3);
            Pindi = Ploadmax(ii)*LoadMix(ii,3)*ppuindi*(1+epsilon*randn);
            Qindi = Pindi*sqrt((1/fpind^2) - 1)*(1+epsilon*randn);
            %Build load matrices
            Pload(ii,kk) = Pcomi + Presi + Pindi;
            Qload(ii,kk) = Qcomi + Qresi + Qindi;
        end
        %Complex aggregated loads
        Sp = diag(Pload(:,kk)+j*Qload(:,kk));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find time in hours
        hr  = floor(t(kk)/60);
        hr  = mod(hr,24);
        if hr>24,
            hr = hr - 24;
        end
        nxthr = min(hr+1,24);
        interhr = mod(t(kk),60)/60;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PV Generation
        Irr = (1-interhr)*PerfilMensualIrr(hr+1) + interhr*PerfilMensualIrr(nxthr+1);
        Psolar(:,kk)  = (unos+epsilon*randn(N,1)).*SPVinstal*Irr/Imax;
        fpsolar(:,kk) = SPfp(kk)*unos;
        Qsolar(:,kk)  = Psolar(:,kk).*sqrt((1./fpsolar(:,kk).^2) - 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RES complex matrix
        SPV = diag(j*Qcomp ... %compensation
                  + Psolar(:,kk) + j*Qsolar(:,kk)); %PV
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
                sss   = SPV*unos - Sp*unos;
                iii   = conj(sss./vvv);
                vvv   = invMtx*iii - bbb ;
                VV    = diag(vvv);

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if norm(h_I(VV,SPV,Sp))<tol,
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
        fp0(kk) = abs( real(s0)/abs(s0) );%*(1+epsilon*randn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New power factor SP for RES
        SPfp(kk+1) = SPfp(kk);
end
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mask for hours with sun
    daytime = sum(Psolar,1)>0;
    mask    = (t>0).*(t<24*60);
    daytime = daytime.*mask;
    tuprise = t(find(daytime,1,'first'))/60;
    tdwrise = t(find(daytime,1,'last'))/60;
    
    k0000 = find(mask,1,'first');
    k2400 = find(mask,1,'last');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % power factor at substation
    figure
        subplot(1,2,2)
            evalc('hold');
            plot(t/60,fp0,'-','LineWidth',1.5)
            title('Open Loop')
            legend({'$f_{p,0}$'},'location','southwest','Interpreter','latex')
            
            xlim([hrs(1),hrs(end)]);
            ylim([fpmin,0.95]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
    % Potencias subestaciOn
        subplot(1,2,1) 
            evalc('hold');
            P0minOL = min(P0(k0000:k2400));
            P0maxOL = max(P0(k0000:k2400));
            plot(t/60,P0*1e-6,'-','LineWidth',1.5)
            title('Open Loop')
            legend({'$P_{0}[MW]$'},'location','northwest','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            ylim([4,10]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
    if grafgrabar,
        nombre   = 'P0_OpenLoop.eps';
        set(gcf,'PaperUnits','centimeters',...
                'PaperSize',[2*ancho alto],...
                'PaperPosition',[0 0 2*ancho alto]); %[0 0 ancho alto]
        print('-depsc','-r200',nombre) % FunciOn para guardar .eps
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % voltages
    figure
            evalc('hold');
            plot(t/60,abs(ve)/V,'-','LineWidth',1.5)
            title('Open Loop')
            legend({'$v_{i}[p.u.]$'},'location','northwest','Interpreter','latex')
            
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potencias carga y PV
    figure
        subplot(1,2,1) 
            evalc('hold');
            plot(t/60,Pload*1e-3,'-','LineWidth',1.5)
            title('$P_{load}[kW]$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(1,2,2) 
            evalc('hold');
            plot(t/60,Psolar*1e-3,'-','LineWidth',1.5)
            title('$P_{PV}[kW]$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
    if grafgrabar,
        nombre   = 'P_Load_PV.eps';
        set(gcf,'PaperUnits','centimeters',...
                'PaperSize',[2*ancho alto],...
                'PaperPosition',[0 0 2*ancho alto]); %[0 0 ancho alto]
        print('-depsc','-r200',nombre) % FunciOn para guardar .eps
    end
    figure
            evalc('hold');
            fpload = Pload./sqrt(Pload.^2 + Qload.^2);
            plot(t/60,fpload,'-','LineWidth',1.5)
            title('$f_{p,load}$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            ylim([0.84,0.90]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hidrogen + consensus
    fprintf('   H2 + Consensus Closed Loop\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consensus + H2 control gains
    fp_H2    = 1;
    
    K11      = 0;%must be zero
    K12      = 0.10;
    K21      = 0.05;
    K22      = 0.1;
    
    K13      = 0;%must be zero
    K23      = 0;
    K31      = 0;%must be zero
    K32      = 0.01;
    K33      = -0.15;%must be negative
    KKK      = [K11, K12, K13;
                K21, K22, K23;
                K31, K32, K33];
    
    fprintf('   Constants KKK= [%1.4f,%1.4f,%1.4f\n                   %1.4f,%1.4f,%1.4f\n                   %1.4f,%1.4f,%1.4f ]\n',...
            K11,K12,K13,...
            K21,K22,K23,...
            K31,K32,K33);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power reference constants
    kinic  = find(t==0,1,'first');
    kfin   = find(t==24*60,1,'first');
    
    E024hrs_1  = sum(P0(kinic:kfin))*h; %W*min
    DPmaxmin   = 2.0e6;
    
    %define rise time of ref
    p0ref_t1 = 17*60;
    p0ref_t2 = 21*60;
    Th       = p0ref_t2 - p0ref_t1;
    
    %define floor and ceil of ref
    p0min    = (E024hrs_1 - DPmaxmin*Th)/(24*60);
    p0max    = p0min + DPmaxmin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some constant matrices
    HHH1 = [0,1,0;
            1,0,0;
            0,0,0]; 
    HHH2 = [0,1;
            0,1;
            1,0]; 
    HHH3 = [0;
            0;
            1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocating vectors
    ve       = zeros(N,length(t));
    vant     = V*unos;

    P0       = zeros(1,length(t));
    Q0       = zeros(1,length(t));
    fp0      = zeros(1,length(t));
    
    fref     = 0.90*ones(1,length(t));
    
    SPfp     = zeros(1,length(t));
    SPfp(1)  = 1;%fpref0;%
    
    P_SP  = zeros(1,length(t));
    PFC_ref  = zeros(1,length(t));
    PH2_ref  = zeros(1,length(t));
    
    PFC      = zeros(N,length(t));
    PH2      = zeros(N,length(t));
    QFC      = zeros(N,length(t));
    QH2      = zeros(N,length(t));
    
    p0ref    = zeros(1,length(t));
    
    MMM      = zeros(2,3,length(t));
    GGG      = zeros(3,3,length(t));
    eigG     = zeros(3,length(t));
    alpha    = zeros(1,length(t));
    epsilon  = zeros(1,length(t));
    
    Sw       = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop on time
for kk=1:length(t),
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ZIP loads already computed
        Sp = diag(Pload(:,kk)+j*Qload(:,kk));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Photo-voltaic Generation (P already computed)
        fpsolar(:,kk) = SPfp(kk)*unos;
        Qsolar(:,kk)  = Psolar(:,kk).*sqrt((1./fpsolar(:,kk).^2) - 1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RES complex matrix
        SPV = diag(j*Qcomp ... %compensation
                  + Psolar(:,kk) + j*Qsolar(:,kk)); %PV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FC injection
        PFC(:,kk) = max(0,min(1,PFC_ref(kk)))*PFCnom;
%         PFC(:,kk) = PFC_ref(kk)*PFCnom;

        fp = SPfp(kk);
        QFC(:,kk) = PFC(:,kk)*sqrt(1/fp^2 - 1);
        
        SFC   = diag(PFC(:,kk) + j*QFC(:,kk));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %H2
        PH2(:,kk) = max(0,min(1,PH2_ref(kk)))*PH2nom;
%         PH2(:,kk) = PH2_ref(kk)*PH2nom;

        fp = fp_H2;
        QH2(:,kk) = PH2(:,kk)*sqrt(1/fp^2 - 1);
        
        SH2   = diag(PH2(:,kk) + j*QH2(:,kk));
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
                sss   = (SPV+SFC)*unos - (Sp+SH2)*unos;
                iii   = conj(sss./vvv);
                vvv   = invMtx*iii - bbb ;
                VV    = diag(vvv);

                cont = cont + 1;
                if cont>=maxit,
                    permiso = 0;
                end
                if norm(h_I(VV,SPV+SFC,Sp+SH2))<tol,
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
        fp0(kk) = abs( real(s0)/abs(s0) );%*(1+epsilon*randn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Built power reference
        if t(kk)<p0ref_t1,
            p0ref(kk) = p0min;
        elseif t(kk)<p0ref_t1+60,
            p0ref(kk) = p0refant + h*(p0max-p0min)/60;
        elseif t(kk)<p0ref_t2,
            p0ref(kk) = p0max;
        elseif t(kk)<p0ref_t2+60,
            p0ref(kk) = p0refant - h*(p0max-p0min)/60;
        else
            p0ref(kk) = p0min;
        end
        p0refant = p0ref(kk);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New set points for RE
        e1    = fref(kk)-fp0(kk);
        e2    = SPfp(kk)-fp0(kk);
        e3    = p0ref(kk) - P0(kk);
        
        temp  = -h*KKK*[e1;e2;e3];
        SPfp(kk+1) = max(fpmin,min(SPfp(kk) +temp(1),1));
        fref(kk+1) = max(fpmin,min(fref(kk) +temp(2),1));
        P_SP(kk+1) = P_SP(kk) +temp(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New H2 and FC references
        HHH0 = [1,0,0;
                0,0,-(1-Sw)/PFCTot;
                0,0, Sw/PH2Tot
                ]; 
        Swant = Sw;
        if P_SP(kk+1)>0
          PFC_ref(kk+1) = 0;
          PH2_ref(kk+1) = P_SP(kk+1)/PH2Tot;
          Sw = 1;
        else
          PFC_ref(kk+1) =-P_SP(kk+1)/PFCTot;
          PH2_ref(kk+1) = 0;
          Sw = 0;
        end              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix M
        aaa = (1-conj(s0)/s0)/(2*abs(s0));
        AAA = V*conj(Y0) + diag(conj(YYY)*conj(vant));
        BBB = diag(vant)*conj(YYY);
        CCC = -j*diag(PFC(:,kk) + Psolar(:,kk))/(SPfp(kk)^2*sqrt(1 - SPfp(kk)^2));
        DDD = (1 + j*sqrt(1/SPfp(kk)^2 - 1))*Id;
        
        MMM0 = 0.5*[1, 1;
                      aaa, conj(aaa)];
        MMM1 = [zeros(1,N), V*unos'*conj(Y0);
                      V*unos'*Y0, zeros(1,N)]/...
                      [AAA, BBB; conj(BBB), conj(AAA)];
        MMM2 = [CCC*unos, DDD*PFCnom, -PH2nom;
                conj(CCC)*unos, conj(DDD)*PFCnom, -PH2nom];
        Maux  = MMM0*MMM1*MMM2;
        MMM(:,:,kk)  = Maux;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matrix G
        Gaux = eye(3,3) - h*(HHH1 - HHH2*Maux*HHH0)*KKK;
        Gaux(isnan(Gaux))=0;
        GGG(:,:,kk)  = Gaux;
        eigG(:,kk)   = eig(Gaux);
        
        alpha(kk)   = Swant*Maux(1,3)/PH2Tot - (1-Swant)*Maux(1,2)/PFCTot;
        epsilon(kk) = Swant*Maux(2,3)/PH2Tot - (1-Swant)*Maux(2,2)/PFCTot;
end
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % power factor at substation
    figure
        subplot(1,2,2)
            evalc('hold');
            plot(t/60,SPfp(1:end-1),'--','LineWidth',1.0)
            plot(t/60,fref(1:end-1),':','LineWidth',1.5)
            plot(t/60,fp0,'-','LineWidth',1.5)
            title('Closed Loop')
            legend({'$f_{p,SP}$','$f_{ref}$','$f_{p,0}$'},'location','southwest','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            ylim([0.85,0.95]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
    % Potencias subestaciOn
        subplot(1,2,1) 
            evalc('hold');
            P0minCL = min(P0(k0000:k2400));
            P0maxCL = max(P0(k0000:k2400));
            plot([t(1) t(end)]/60,[1 1]*p0min*1e-6,':k','LineWidth',1.0)
            plot([t(1) t(end)]/60,[1 1]*p0max*1e-6,':k','LineWidth',1.0)
            pP0ref = plot(t/60,p0ref*1e-6,':r','LineWidth',1.5);
            pP0    = plot(t/60,P0*1e-6,'-','LineWidth',1.5);
            title('Closed loop')
            legend([pP0,pP0ref],{'$P_{0}[MW]$','$P_{ref}[MW]$'},'location','northwest','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            ylim([4,10]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
    if grafgrabar,
        nombre   = 'P0_ClosedLoop.eps';
        set(gcf,'PaperUnits','centimeters',...
                'PaperSize',[2*ancho alto],...
                'PaperPosition',[0 0 2*ancho alto]); %[0 0 ancho alto]
        print('-depsc','-r200',nombre) % FunciOn para guardar .eps
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % voltages
    figure
            evalc('hold');
            plot(t/60,abs(ve)/V,'-','LineWidth',1.5)
            title('Closed Loop')
            legend({'$v_{i}[p.u.]$'},'location','northwest','Interpreter','latex')
            
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% approximated and actual power at substation
    figure
        s0approx = unos'*((Pload+j*Qload + PH2+j*QH2) - (PFC+j*QFC + Psolar+j*Qsolar));
        subplot(2,2,1)
            evalc('hold');
            plot(t/60,P0*1e-6,'-','LineWidth',1.5);
            plot(t/60,real(s0approx)*1e-6,'-','LineWidth',1.5);
            title('Closed loop')
            legend({'$P_{0}[MW]$','$P_{0,approx}[MW]$'},'location','northwest','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,2,2)
            evalc('hold');
            plot(t/60,Q0*1e-6,'-','LineWidth',1.5);
            plot(t/60,imag(s0approx)*1e-6,'-','LineWidth',1.5);
            title('Closed loop')
            legend({'$Q_{0}[MW]$','$Q_{0,approx}[MW]$'},'location','northwest','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,2,3)
            evalc('hold');
            plot(t/60,100*(P0-real(s0approx))./P0,'-','LineWidth',1.5);
            legend({'$100(P_0 - P_{0,approx})/P_{0}$'},'location','northwest','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,2,4)
            evalc('hold');
            plot(t/60,100*(Q0-imag(s0approx))./Q0,'-','LineWidth',1.5);
            legend({'$100(Q_0 - Q_{0,approx})/Q_{0}$'},'location','northwest','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FC-H2 power
    figure
        subplot(1,2,1) 
            evalc('hold');
            plot(t/60,PH2*1e-6,'-','LineWidth',1.5)
            title('$P_{H2}[MW]$','Interpreter','latex')
%             legend({'$P_{H2}[MW]$'},'location','northeast','Interpreter','latex'))
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);

        subplot(1,2,2) 
            evalc('hold');
            plot(t/60,PFC*1e-6,'-','LineWidth',1.5)
            title('$P_{FC}[MW]$','Interpreter','latex')
%             legend({'$P_{FC}[MW]$'},'location','northeast','Interpreter','latex'))
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
    if grafgrabar,
        nombre   = 'H2FC_ClosedLoop.eps';
        set(gcf,'PaperUnits','centimeters',...
                'PaperSize',[2*ancho alto],...
                'PaperPosition',[0 0 2*ancho alto]); %[0 0 ancho alto]
        print('-depsc','-r200',nombre) % FunciOn para guardar .eps
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % matrix MMM
    figure
        subplot(2,3,1) 
            evalc('hold');
            plot(t/60,squeeze(real(MMM(1,1,:))),'-','LineWidth',1.5)
            title('$[M]_{11}$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,3,2) 
            evalc('hold');
            plot(t/60,ones(size(t))*-PFCTot,':','LineWidth',1.5)
            plot(t/60,squeeze(real(MMM(1,2,:))),'-','LineWidth',1.5)
            title('$[M]_{12}$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,3,3) 
            evalc('hold');
            plot(t/60,ones(size(t))*PH2Tot,':','LineWidth',1.5)
            plot(t/60,squeeze(real(MMM(1,3,:))),'-','LineWidth',1.5)
            title('$[M]_{13}$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,3,4) 
            evalc('hold');
            hatM21 = -sum(Psolar+PFC,1)./P0;
            hatM21min = -PPVmax/p0min;
            plot(t/60,hatM21,':','LineWidth',1.5)
            plot(t/60,squeeze(real(MMM(2,1,:))),'-','LineWidth',1.5)
            plot(t/60,hatM21min*ones(size(t)),'k:','LineWidth',1.0)
            title('$[M]_{21}$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            ylim([-0.8,0]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,3,5) 
            evalc('hold');
            hatM22 = 0.5*fp0.*(1-fp0.^2)*PFCTot./P0;
            plot(t/60,hatM22,':','LineWidth',1.5)
            plot(t/60,squeeze(real(MMM(2,2,:))),'-','LineWidth',1.5)
            title('$[M]_{22}$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(2,3,6) 
            evalc('hold');
            hatM23 = 0.5*fp0.*(1-fp0.^2)*PH2Tot./P0;
            plot(t/60,hatM23,':','LineWidth',1.5)
            plot(t/60,squeeze(real(MMM(2,3,:))),'-','LineWidth',1.5)
            title('$[M]_{23}$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % matrix GGG
    Fontsize = 12;
    figure
        subplot(3,3,1) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(1,1,:))),'-','LineWidth',1.5)
            plot(t/60,real(eigG(1,:)),':','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{11}$','$\lambda_1$'},...
                    'location','south','Orientation','horizontal','Interpreter','latex')
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,2) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(1,2,:))),'-','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{12}$'},...
                    'location','north','Orientation','horizontal','Interpreter','latex')
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,3) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(1,3,:))),'-','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{13}$'},...
                    'location','south','Orientation','horizontal','Interpreter','latex')
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,4) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(2,1,:))),'-','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{21}$'},...
                    'location','south','Orientation','horizontal','Interpreter','latex')
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,5) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(2,2,:))),'b-','LineWidth',1.5)
            plot(t/60,real(eigG(2,:)),'r--','LineWidth',1.5)
            plot(t/60,1+h*K12*(hatM21-1),'g:','LineWidth',1.5)
            plot(t/60,(1+h*K12*(hatM21min-1))*ones(size(t)),'k:','LineWidth',1.0)
            plot(t/60,(1+h*K12*(0-1))*ones(size(t)),'k:','LineWidth',1.0)
            xlim([hrs(1),hrs(end)]);
            ylim([0.2,0.7]);
            legend({'$G_{22}$','$\lambda_2$','$1+hK_{12}\hat{M}_{21}$'},...
                    'location','north','Orientation','vertical','Interpreter','latex')
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,6) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(2,3,:))),'-','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{23}$'},...
                    'location','south','Orientation','horizontal','Interpreter','latex')
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,7) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(3,1,:))),'-','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{31}$'},...
                    'location','south','Orientation','horizontal','Interpreter','latex')
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,8) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(3,2,:))),'-','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{32}$'},...
                    'location','south','Orientation','horizontal','Interpreter','latex')
            box ON
            grid ON
%             ylim([-3e4,1e4])
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(3,3,9) 
            evalc('hold');
            plot(t/60,squeeze(real(GGG(3,3,:))),'-b','LineWidth',1.5)
            plot(t/60,real(eigG(3,:)),'--r','LineWidth',1.5)
            plot([t(1),t(end)]/60,[1 1]*(1+h*K33),'g:','LineWidth',1.5)
            xlim([hrs(1),hrs(end)]);
            legend({'$G_{33}$','$\lambda_3$','$1+hK_{33}$'},...
                    'location','north','Orientation','vertical','Interpreter','latex')
            box ON
            grid ON
            ylim([0.22,0.29]);
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtxtsh)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize*0.75);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
    if grafgrabar,
        nombre   = 'Gmatrix_ClosedLoop.eps';
        set(gcf,'PaperUnits','centimeters',...
                'PaperSize',[3*ancho 3*alto],...
                'PaperPosition',[0 0 3*ancho 3*alto]); %[0 0 ancho alto]
        print('-depsc','-r200',nombre) % FunciOn para guardar .eps
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % alpha and epsilon
    figure
        subplot(1,2,1) 
            evalc('hold');
            plot(t/60,abs(alpha),'-','LineWidth',1.5)
            title('$\alpha$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
        subplot(1,2,2) 
            evalc('hold');
            plot(t/60,abs(epsilon),'-','LineWidth',1.5)
            title('$\epsilon$','Interpreter','latex')
            xlim([hrs(1),hrs(end)]);
            box ON
            grid ON
            hAx=gca;  % avoid repetitive function calls
                set(hAx,'XTick',hrs)
                set(hAx,'XTickLabel',hrtext)
                set(hAx,'xminorgrid','off','yminorgrid','off')
                set(hAx,'Layer','top');
                set(gca,'FontName','Times New Roman','FontSize',Fontsize);
                ylimits = ylim;
                rect=rectangle('Position',[tuprise,ylimits(1),tdwrise-tuprise,ylimits(2)-ylimits(1)],...
                               'FaceColor',[255 255 230]/255,...
                               'EdgeColor',[20 20 20]/255,'LineStyle',':');
                uistack(rect,'bottom') ;
                ylim(ylimits);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('   Power constant values\n')
        fprintf('    Instaled C Shunt: %f [kVAr]\n',sum(Qcomp)*1e-3)
        fprintf('    Instaled PV gen: %f [MW]\n',PPVmax*1e-6)
        fprintf('    Instaled Electrolizers: %f [MW]\n',PH2Tot*1e-6)
        fprintf('    Instaled FC gen: %f [MW]\n',PFCTot*1e-6)
        fprintf('    Max. PV gen: %f [MW]\n',max(sum(Psolar,1)*1e-6))
        fprintf('    Max. P load: %f [MW]\n',max(sum(Pload,1)*1e-6))
        fprintf('    Min. P load: %f [MW]\n',min(sum(Pload,1)*1e-6))
        fprintf('    Max.PV gen / Min. P load ratio: %f [%%]\n',100*max(sum(Psolar,1))/min(sum(Pload,1)))
        fprintf('    Open Loop Min. P0: %f [MW]\n',P0minOL*1e-6)
        fprintf('    Open Loop Max. P0: %f [MW]\n',P0maxOL*1e-6)
        fprintf('    Min. P0 ref: %f [MW]\n',p0min*1e-6)
        fprintf('    Max. P0 ref: %f [MW]\n',p0max*1e-6)
        fprintf('    Closed Loop Min. P0: %f [MW]\n',P0minCL*1e-6)
        fprintf('    Closed Loop Max. P0: %f [MW]\n',P0maxCL*1e-6)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('   Slope approximations\n')
        M21min = -PPVmax/p0min;
        M21max = 0;
        fprintf('    Installed PV gen. - Min. P0 ref ratio: %f\n',-M21min)
        fprintf('    Lower bound for G22: %f\n',1+h*K12*(M21min-1))
        fprintf('    Upper bound for G22: %f\n',1+h*K12*(M21max-1))
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('   Comparison between cases\n')
    E024hrs_2  = sum(P0(kinic:kfin))*h; %W*min
        fprintf('    Daily energy injection at substation:\n')
        fprintf('      Open loop: %f [MWhr]\n',E024hrs_1*1e-6/60)
        fprintf('      Closed loop: %f [MWhr]\n',E024hrs_2*1e-6/60)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EH2    = h*sum(sum(PH2(:,kinic:kfin),1))/60;
    EFC    = h*sum(sum(PFC(:,kinic:kfin),1))/60;
    EPV    = h*sum(sum(Psolar(:,kinic:kfin),1))/60;
    Eload  = h*sum(sum(Pload(:,kinic:kfin),1))/60;
        fprintf('    Daily energy converted to hydrogen: %f [MWhr]\n',EH2*1e-6)
        fprintf('    Daily energy injected by fuel cells: %f [MWhr]\n',EFC*1e-6)
        fprintf('    Net hydrogen gain: %f [kWhr]\n',(EH2-EFC)*1e-3)
        fprintf('    Daily energy injected by PV: %f [MWhr]\n',EPV*1e-6)
        fprintf('    Load daily energy demand: %f [MWhr]\n',Eload*1e-6)
        fprintf('    daily energy PV/Load ratio: %f [%%]\n',EPV/Eload*100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
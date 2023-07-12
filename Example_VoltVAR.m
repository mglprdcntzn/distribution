%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for distribution with DG Volt-VAR + C
% Miguel Parada Contzen
% Conce, julio 2023
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
    mes    = 1; %month of the year
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
% VoltVAR
    fprintf('   Volt-Var\n')
    K_VoltVAR = 10;
    p         = 3/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocating vectors
    Pload  = zeros(N,length(t));
    Qload  = zeros(N,length(t));

    Psolar   = zeros(N,length(t));
    fpsolar  = ones(N,length(t));
    Qsolar   = zeros(N,length(t));

    ve       = zeros(N,length(t));
    vant     = V*unos;
    vemean   = vant;

    P0       = zeros(1,length(t));
    Q0       = zeros(1,length(t));
    fp0      = zeros(1,length(t));
    
    SPfp     = zeros(N,length(t));
    SPfp(:,1)  = ones(N,1);%fpref0;%    
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
        fpsolar(:,kk) = (Psolar(:,kk)>0).*SPfp(:,kk) + (Psolar(:,kk)<=0);
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
        vemean   = p*vemean + (1-p)*abs(vant); %mean voltage approx
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find powers at substation root node 0
        s0      = conj(Y00)*V^2 + V*unos'*conj(Y0)*conj(vant);
        P0(kk)  = real(s0);
        Q0(kk)  = imag(s0);
        fp0(kk) = abs( real(s0)/abs(s0) );%*(1+epsilon*randn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New power factor SP for RES
        eve = (abs(vant) - vemean)./vemean; %measured voltages error w/r to mean
        SPfp(:,kk+1) = max(fpmin,min(1,SPfp(:,kk) - K_VoltVAR*eve));
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
            title('Volt-VAR')
            legend({'$f_{p,0}$'},'location','southwest','Interpreter','latex')
            
            xlim([hrs(1),hrs(end)]);
            ylim([fpmin,1]);
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
            title('Volt-VAR')
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % voltages
    figure
            evalc('hold');
            plot(t/60,abs(ve)/V,'-','LineWidth',1.5)
            title('Volt-VAR')
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
    % Potencias carga
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
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potencias PV
    figure
        subplot(1,2,1) 
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
        subplot(1,2,2) 
            evalc('hold');
            plot(t/60,fpsolar,'-','LineWidth',1.5)
            title('$fp_{PV}$','Interpreter','latex')
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
    % power factor at substation
    figure
        subplot(1,2,1)
            evalc('hold');
            plot(t/60,fp0,'-','LineWidth',1.5)
            title('Volt-VAR: $f_{p,0}$','Interpreter','latex')
            %legend({'$f_{p,0}$'},'location','southwest','Interpreter','latex')
            
            xlim([hrs(1),hrs(end)]);
            ylim([fpmin,1]);
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
                ylim(ylimits);subplot(1,2,2) 
        subplot(1,2,2)
            evalc('hold');
            plot(t/60,fpsolar,'-','LineWidth',1.5)
            title('Volt-VAR: $f_{p,PV}$','Interpreter','latex')
            %title('$fp_{PV}$','Interpreter','latex')
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
        nombre   = 'fp_VoltVAR.eps';
        set(gcf,'PaperUnits','centimeters',...
                'PaperSize',[2*ancho alto],...
                'PaperPosition',[0 0 2*ancho alto]); %[0 0 ancho alto]
        print('-depsc','-r200',nombre) % FunciOn para guardar .eps
    end
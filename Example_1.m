%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1 for distribution
% Miguel Parada Contzen
% Conce, diciembre 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
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
%Static arbitrary load 
    fp    = 0.90;%equal pwr fctr for everyone
    Ppu   = 0.4; % per-unit of max load
    Sload = Ppu*diag(NodeData(2:end,4))*1e3*( 1 + j*sqrt(1/fp^2 - 1));
%Static RES
    Sres  = 0.25*diag(NodeData(2:end,8)) +...%PV
            0.30*diag(NodeData(2:end,9)) +...%Wind
               j*diag(NodeData(2:end,10));%compensation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Branches
    BranchData = dlmread('branches.csv',';',1,0);
    NB = size(BranchData,1);
    
    DDD = zeros(NB,N+1); %incidence matrix
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
% Algorithm tunning
tol        = 1e-10; %tolerance
maxit      = 1000; %max number of iteration
terminator = 0;%0: f_I; 1: f_P
    fprintf('   N° iteracione maximas: %i\n',maxit)
    fprintf('   Tolerancia: %e\n',tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exponential function (expm too slow?)
my_expmj = @(XXX) diag(  cos(diag(XXX)) + j*sin(diag(XXX)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power balance functions
h_P = @(VVV)   Sres*unos - Sload*unos - V*conj(Y0)*VVV*unos - VVV*conj(YYY)*conj(VVV)*unos;
A_P = @(VVV)  -V*conj(Y0) - diag(conj(YYY)*conj(VVV)*unos);
B_P = @(VVV)  -VVV*conj(YYY);
C_P = @(AAA,BBB,RRR,Psi)  AAA*my_expmj(Psi) + BBB*my_expmj(-Psi);
D_P = @(AAA,BBB,RRR,Psi)  AAA*j*RRR*my_expmj(Psi) - BBB*j*RRR*my_expmj(-Psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current balance functions
h_I = @(VVV)  VVV\Sres*unos - VVV\Sload*unos - conj(Y0)*unos*V - conj(YYY)*conj(VVV)*unos;
A_I = @(VVV)  -(VVV^2)\Sres + (VVV^2)\Sload;
B_I = @(VVV)  -conj(YYY);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terminator functions
if terminator,
    TheEnd = @(VVV) norm(h_P(VVV))<tol;
    fprintf('   Terminator: h_P\n')
else
    TheEnd = @(VVV) norm(h_I(VVV))<tol;
    fprintf('   Terminator: h_I\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation initial conditions
V0   = V*eye(N,N);
R0   = V0;
PSI0 = zeros(N,N);
uw0 = [real(V0)*unos; imag(V0)*unos];
rp0 = [R0*unos; PSI0*unos];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITMO 1: NR-clAsico
fprintf('ALGORITMO 1: NR-clAsico\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condiciones iniciales
    VV   = V0;
    RR   = R0;
    PSI  = PSI0;
    rp   = rp0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iteraciones algoritmo NR-clAsico
    permiso = 1;
    cont    = 0;
    tic;
    while permiso,
        h   = h_P(VV);

        AA  = A_P(VV);
        BB  = B_P(VV);
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
        if TheEnd(VV),
            permiso = 0;
        end
    end
    elapsedTime = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guardar voltajes resultantes
    ve      = VV*unos;
    veAmpu  = abs(ve)/V;
    vePhdeg = angle(ve)*180/pi;
    ve_1    = ve;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % presentar resultados
    figure(1)
    subplot(2,3,1)
    plot([1;veAmpu],'*-')
    grid ON
    
    figure(2)
    subplot(2,3,1)
    plot([0;vePhdeg],'*-')
    grid ON
    
    fprintf('   N° iteraciones: %i\n',cont)
    fprintf('   Tiempo: %1.1f [ms]\n',elapsedTime*1000)
    fprintf('   ||h_I||: %1.4e\n',norm(h_I(VV)))
    output = sprintf('& %1.1f & %i & %1.2e \n', elapsedTime*1000, cont, norm(h_I(VV)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITMO 2: NR-UW
fprintf('ALGORITMO 2: NR-UW\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condiciones iniciales
    VV   = V0;
    uw   = uw0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iteraciones algoritmo NR-UW
    permiso = 1;
    cont    = 0;
    tic;
    while permiso,
        h   = h_P(VV);

        AA  = A_P(VV);
        BB  = B_P(VV);

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
        if TheEnd(VV),
            permiso = 0;
        end
    end
    elapsedTime = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guardar voltajes resultantes
    ve      = VV*unos;
    veAmpu  = abs(ve)/V;
    vePhdeg = angle(ve)*180/pi;
    ve_2    = ve;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % presentar resultados
    figure(1)
    subplot(2,3,2)
    plot([1;veAmpu],'*-')
    grid ON
    
    figure(2)
    subplot(2,3,2)
    plot([0;vePhdeg],'*-')
    grid ON
    
    fprintf('   N° iteraciones: %i\n',cont)
    fprintf('   Tiempo: %1.1f [ms]\n',elapsedTime*1000)
    fprintf('   ||h_I||: %1.4e\n',norm(h_I(VV)))
    output = [output, sprintf('& %1.1f & %i & %1.2e \n', elapsedTime*1000, cont, norm(h_I(VV)))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITMO 3: NR-IUW
fprintf('ALGORITMO 3: NR-IUW\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condiciones iniciales
    VV   = V0;
    uw   = uw0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iteraciones algoritmo NR-IUW
    permiso = 1;
    cont    = 0;
    tic;
    while permiso,
        h   = h_I(VV);

        AA  = A_I(VV);
        BB  = B_I(VV);

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
        if TheEnd(VV),
            permiso = 0;
        end
    end
    elapsedTime = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guardar voltajes resultantes
    ve      = VV*unos;
    veAmpu  = abs(ve)/V;
    vePhdeg = angle(ve)*180/pi;
    ve_3    = ve;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % presentar resultados
    figure(1)
    subplot(2,3,3)
    plot([1;veAmpu],'*-')
    grid ON
    
    figure(2)
    subplot(2,3,3)
    plot([0;vePhdeg],'*-')
    grid ON
    
    fprintf('   N° iteraciones: %i\n',cont)
    fprintf('   Tiempo: %1.1f [ms]\n',elapsedTime*1000)
    fprintf('   ||h_I||: %1.4e\n',norm(h_I(VV)))
    output = [output, sprintf('& %1.1f & %i & %1.2e \n', elapsedTime*1000, cont, norm(h_I(VV)))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITMO 4: NR-IUWsimpl
fprintf('ALGORITMO 4: NR-IUWsimpl\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condiciones iniciales
    VV   = V0;
    uw   = uw0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iteraciones algoritmo NR-UW
    permiso = 1;
    cont    = 0;
    tic;
    Mtx = [-real(YYY), imag(YYY);
            imag(YYY), real(YYY)];
    invMtx = Mtx\eye(size(Mtx));
    elapsedTimeMtx = toc;
    tic;
    while permiso,
        h   = h_I(VV);
        
        vtr = -[real(h);
                imag(h)];

        Delta = invMtx*vtr;
        uw    = uw + Delta;
        VV    = diag(uw(1:N)) + j*diag(uw(N+1:end));

        cont = cont + 1;
        if cont>=maxit,
            permiso = 0;
        end
        if TheEnd(VV),
            permiso = 0;
        end
    end
    elapsedTime = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guardar voltajes resultantes
    ve      = VV*unos;
    veAmpu  = abs(ve)/V;
    vePhdeg = angle(ve)*180/pi;
    ve_4    = ve;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % presentar resultados
    figure(1)
    subplot(2,3,4)
    plot([1;veAmpu],'*-')
    grid ON
    
    figure(2)
    subplot(2,3,4)
    plot([0;vePhdeg],'*-')
    grid ON
    
    fprintf('   N° iteraciones: %i\n',cont)
    fprintf('   Tiempo Mtx: %1.1f [ms]\n',elapsedTimeMtx*1000)
    fprintf('   Tiempo iter: %1.1f [ms]\n',elapsedTime*1000)
    fprintf('   Tiempo total: %1.1f [ms]\n',(elapsedTime+elapsedTimeMtx)*1000)
    fprintf('   ||h_I||: %1.4e\n',norm(h_I(VV)))
    output = [output, sprintf('& %1.1f & %i & %1.2e \n', elapsedTime*1000, cont, norm(h_I(VV)))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITMO 5: NR-Corrientes sucesivas
fprintf('ALGORITMO 5: NR-Corrientes sucesivas\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % condiciones iniciales
    VV   = V0;
    ve   = VV*unos;
    sss  = (Sres-Sload)*unos;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iteraciones algoritmo NR-Corrientes sucesivas
    permiso = 1;
    cont    = 0;
    tic;
    invMtx = YYY\eye(size(YYY));
    bbb    = invMtx*Y0*unos*V;
    elapsedTimeMtx = toc;
    tic;
    while permiso,
        iii   = conj(sss./ve);
        ve   = invMtx *iii - bbb ;
        VV    = diag(ve);

        cont = cont + 1;
        if cont>=maxit,
            permiso = 0;
        end
        if TheEnd(VV),
            permiso = 0;
        end
    end
    elapsedTime = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Guardar voltajes resultantes
    ve      = VV*unos;
    veAmpu  = abs(ve)/V;
    vePhdeg = angle(ve)*180/pi;
    ve_5    = ve;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % presentar resultados segUn numeraciOn paper
    figure(1)
    subplot(2,3,5)
    plot([1;veAmpu],'*-')
    grid ON
    
    figure(2)
    subplot(2,3,5)
    plot([0;vePhdeg],'*-')
    grid ON
    
    fprintf('   N° iteraciones: %i\n',cont)
    fprintf('   Tiempo Mtx: %1.1f [ms]\n',elapsedTimeMtx*1000)
    fprintf('   Tiempo iter: %1.1f [ms]\n',elapsedTime*1000)
    fprintf('   Tiempo total: %1.1f [ms]\n',(elapsedTime+elapsedTimeMtx)*1000)
    fprintf('   ||h_I||: %1.4e\n',norm(h_I(VV)))
    output = [output, sprintf('& %1.1f & %i & %1.2e \n', elapsedTime*1000, cont, norm(h_I(VV)))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%presentar todos los voltajes
    veAmpu_1  = abs(ve_1)/V;
    vePhdeg_1 = angle(ve_1)*180/pi;
    veAmpu_2  = abs(ve_2)/V;
    vePhdeg_2 = angle(ve_2)*180/pi;
    veAmpu_3  = abs(ve_3)/V;
    vePhdeg_3 = angle(ve_3)*180/pi;
    veAmpu_4  = abs(ve_4)/V;
    vePhdeg_4 = angle(ve_4)*180/pi;
    veAmpu_5  = abs(ve_5)/V;
    vePhdeg_5 = angle(ve_5)*180/pi;
    
    e2Am_12 = ([1;veAmpu_1]-[1;veAmpu_2]).^2;
    fprintf('   Error Am acumulado 1-2: %e\n',sum(e2Am_12))
    e2Am_13 = ([1;veAmpu_1]-[1;veAmpu_3]).^2;
    fprintf('   Error Am acumulado 1-3: %e\n',sum(e2Am_13))
    e2Am_14 = ([1;veAmpu_1]-[1;veAmpu_4]).^2;
    fprintf('   Error Am acumulado 1-4: %e\n',sum(e2Am_14))
    e2Am_15 = ([1;veAmpu_1]-[1;veAmpu_5]).^2;
    fprintf('   Error Am acumulado 1-5: %e\n',sum(e2Am_15))
    
    e2Ph_12 = ([1;vePhdeg_1]-[1;vePhdeg_2]).^2;
    e2Ph_13 = ([1;vePhdeg_1]-[1;vePhdeg_3]).^2;
    e2Ph_14 = ([1;vePhdeg_1]-[1;vePhdeg_4]).^2;
    e2Ph_15 = ([1;vePhdeg_1]-[1;vePhdeg_5]).^2;
    
    ii = 1:N+1;
    % presentar resultados de todos segUn numeraciOn paper
    figure(1)
    subplot(2,3,6)
    plot(ii,e2Am_12,'-*',ii,e2Am_13,'-*',ii,e2Am_14,'-*',ii,e2Am_15,'-*')
    grid ON
    
    figure(2)
    subplot(2,3,6)
    plot(ii,e2Ph_12,'-*',ii,e2Ph_13,'-*',ii,e2Ph_14,'-*',ii,e2Ph_15,'-*')
    grid ON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GrAfico de voltajes segUn NR clAsico
    ancho = 200;
    alto  = 200;
    Fontsize = 12;
    
    figure
    subplot(1,2,1)
        plot(0:N,[1;abs(ve_1)/V],'-*')
            box ON
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            title('|v| [p.u.]')
            xlabel('bus N°');
            xlim([0,N]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    subplot(1,2,2)
        plot(0:N,[0;angle(ve_1)*180/pi],'-*')
            box ON
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            title('\angle v [°]')
            xlabel('bus N°');
            xlim([0,N]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
            
    set(gcf,'PaperUnits','centimeters',...
            'PaperSize',[2*ancho 1*alto],...
            'PaperPosition',[0 0 2*ancho 1*alto]); %[0 0 ancho alto]
    print('-depsc','-r200','v_NRC') % FunciOn para guardar .eps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%imprimir output
disp(output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

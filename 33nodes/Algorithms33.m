%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIP-I2: Algoritmos para resolver ctos de distribuciOn
% Miguel Parada Contzen
% Conce, mayo 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINICIONES PARA P1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definiciones grales
j = sqrt(-1);
V = 12.66e3;
N = 33-1; %N of buses
NB = 35-3; %N of branches

unos = ones(N,1);
Id   = eye(N,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TopologIa de circuito
DDD = zeros(NB,N+1);
DDD(01,01) = 1; DDD(01,02) = -1; %branch 1: 1- 2
DDD(02,02) = 1; DDD(02,03) = -1; %branch 2: 2- 3
DDD(03,03) = 1; DDD(03,04) = -1; %branch 3: 3- 4
DDD(04,04) = 1; DDD(04,05) = -1; %branch 4: 4- 5
DDD(05,05) = 1; DDD(05,06) = -1; %branch 5: 5- 6
DDD(06,06) = 1; DDD(06,07) = -1; %branch 6: 6- 7
DDD(07,07) = 1; DDD(07,08) = -1; %branch 7: 7- 8
DDD(08,08) = 1; DDD(08,09) = -1; %branch 8: 8- 9
DDD(09,09) = 1; DDD(09,10) = -1; %branch 9: 9- 10

DDD(10,10) = 1; DDD(10,11) = -1; %branch 10: 10- 11
DDD(11,11) = 1; DDD(11,12) = -1; %branch 11: 11- 12
DDD(12,12) = 1; DDD(12,13) = -1; %branch 12: 12- 13
DDD(13,13) = 1; DDD(13,14) = -1; %branch 13: 13- 14
DDD(14,14) = 1; DDD(14,15) = -1; %branch 14: 14- 15
DDD(15,15) = 1; DDD(15,16) = -1; %branch 15: 15- 16
DDD(16,16) = 1; DDD(16,17) = -1; %branch 16: 16- 17
DDD(17,17) = 1; DDD(17,18) = -1; %branch 17: 17- 18
DDD(18,02) = 1; DDD(18,19) = -1; %branch 18:  2- 19
DDD(19,19) = 1; DDD(19,20) = -1; %branch 19: 19- 20

DDD(20,20) = 1; DDD(20,21) = -1; %branch 20: 20- 21
DDD(21,21) = 1; DDD(21,22) = -1; %branch 21: 21- 22
DDD(22,03) = 1; DDD(22,23) = -1; %branch 22:  3- 23
DDD(23,23) = 1; DDD(23,24) = -1; %branch 23: 23- 24
DDD(24,24) = 1; DDD(24,25) = -1; %branch 24: 24- 25
DDD(25,06) = 1; DDD(25,26) = -1; %branch 25:  6- 16
DDD(26,26) = 1; DDD(26,27) = -1; %branch 26: 26- 27
DDD(27,27) = 1; DDD(27,28) = -1; %branch 27: 27- 28
DDD(28,28) = 1; DDD(28,29) = -1; %branch 28: 28- 29
DDD(29,29) = 1; DDD(29,30) = -1; %branch 29: 29- 30

DDD(30,30) = 1; DDD(30,31) = -1; %branch 30: 30- 31
DDD(31,31) = 1; DDD(31,32) = -1; %branch 31: 31- 32
DDD(32,32) = 1; DDD(32,33) = -1; %branch 32: 32- 33

% DDD(33,21) = 1; DDD(33,08) = -1; %sw branch 33: 21-  8
% DDD(34,12) = 1; DDD(34,22) = -1; %sw branch 34: 12- 22
% DDD(35,25) = 1; DDD(35,29) = -1; %sw branch 35: 25- 29
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%impendancias de lInea
WWW = diag(1./[
        0.0922 + 0.0470*j; % branch 1
        0.4930 + 0.2511*j; % branch 2
        0.3660 + 0.1864*j; % branch 3
        0.3811 + 0.1941*j; % branch 4
        0.8190 + 0.7070*j; % branch 5
        0.1872 + 0.6188*j; % branch 6
        0.7114 + 0.2351*j; % branch 7
        1.0300 + 0.7400*j; % branch 8
        1.0440 + 0.7400*j; % branch 9
        0.1966 + 0.0650*j; % branch 10
        0.3744 + 0.1238*j; % branch 11
        1.4680 + 1.1550*j; % branch 12
        0.5416 + 0.7129*j; % branch 13
        0.5910 + 0.5260*j; % branch 14
        0.7463 + 0.5450*j; % branch 15
        1.2890 + 1.7210*j; % branch 16
        0.7320 + 0.5740*j; % branch 17
        0.1640 + 0.1565*j; % branch 18
        1.5042 + 1.3554*j; % branch 19
        0.4095 + 0.4784*j; % branch 20
        0.7089 + 0.9373*j; % branch 21
        0.4512 + 0.3083*j; % branch 22
        0.8980 + 0.7091*j; % branch 23
        0.8960 + 0.7011*j; % branch 24
        0.2030 + 0.1034*j; % branch 25
        0.2842 + 0.1447*j; % branch 26
        1.0590 + 0.9337*j; % branch 27
        0.8042 + 0.7006*j; % branch 28
        0.5075 + 0.2585*j; % branch 29
        0.9744 + 0.9630*j; % branch 30
        0.3105 + 0.3619*j; % branch 31
        0.3410 + 0.5302*j%; % branch 32
%         2.0000 + 2.0000*j; % sw branch 33
%         2.0000 + 2.0000*j; % sw branch 34
%         0.5000 + 0.5000*j; % sw branch 35
    ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrices 
hatYYY = DDD'*WWW*DDD;
Y00    = hatYYY(1,1);
Y0     = diag(hatYYY(1,2:N+1));
YYY    = hatYYY(2:end,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Potencias nominales de carga
Sload = diag([
        0.100 + 0.060*j; %node 1
        0.090 + 0.040*j; %node 2
        0.120 + 0.080*j; %node 3
        0.060 + 0.030*j; %node 4
        0.060 + 0.020*j; %node 5
        0.200 + 0.100*j; %node 6
        0.200 + 0.100*j; %node 7
        0.060 + 0.020*j; %node 8
        0.060 + 0.020*j; %node 9
        0.045 + 0.030*j; %node 10
        0.060 + 0.035*j; %node 11
        0.060 + 0.035*j; %node 12
        0.120 + 0.080*j; %node 13
        0.060 + 0.010*j; %node 14
        0.060 + 0.020*j; %node 15
        0.060 + 0.020*j; %node 16
        0.090 + 0.040*j; %node 17
        0.090 + 0.040*j; %node 18
        0.090 + 0.040*j; %node 19
        0.090 + 0.040*j; %node 20
        0.090 + 0.040*j; %node 21
        0.090 + 0.050*j; %node 22
        0.420 + 0.200*j; %node 23
        0.420 + 0.200*j; %node 24
        0.060 + 0.025*j; %node 25
        0.060 + 0.025*j; %node 26
        0.060 + 0.020*j; %node 27
        0.120 + 0.070*j; %node 28
        0.200 + 0.600*j; %node 29
        0.150 + 0.070*j; %node 30
        0.210 + 0.100*j; %node 31
        0.060 + 0.040*j; %node 32
        ])*1e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GeneraciOn RES + reactive compensation
Sres = zeros(N,N);
node = 18-1;
Sres(node,node) = 0.2e6 + j*0.4e6;
node = 22-1;
Sres(node,node) = 0.2e6 ;
node = 25-1;
Sres(node,node) = 0.2e6 ;
node = 33-1;
Sres(node,node) = 0.2e6 + j*0.6e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINICIONES PARA P2 y P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SintonizaciOn algoritmos
tol        = 1e-5; %tolerancia de iteraciones
maxit      = 1000; %N max iteraciones
terminator = 0;%0: h_I; 1: h_P
    fprintf('   No. iteraciones maximas: %i\n',maxit)
    fprintf('   Tolerancia: %e\n',tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FunciOn exponencial (expm es muy lento?)
my_expmj = @(XXX) diag(  cos(diag(XXX)) + j*sin(diag(XXX)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funciones balance de Potencia
h_P = @(VVV)   Sres*unos - Sload*unos - V*conj(Y0)*VVV*unos - VVV*conj(YYY)*conj(VVV)*unos;
A_P = @(VVV)  -V*conj(Y0) - diag(conj(YYY)*conj(VVV)*unos);
B_P = @(VVV)  -VVV*conj(YYY);
C_P = @(AAA,BBB,RRR,Psi)  AAA*my_expmj(Psi) + BBB*my_expmj(-Psi);
D_P = @(AAA,BBB,RRR,Psi)  AAA*j*RRR*my_expmj(Psi) - BBB*j*RRR*my_expmj(-Psi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funciones balance de Corriente
h_I = @(VVV)  VVV\Sres*unos - VVV\Sload*unos - conj(Y0)*unos*V - conj(YYY)*conj(VVV)*unos;
A_I = @(VVV)  -(VVV^2)\Sres + (VVV^2)\Sload;
B_I = @(VVV)  -conj(YYY);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funciones terminator
if terminator,
    TheEnd = @(VVV) norm(h_P(VVV))<tol;
    fprintf('   Terminator: h_P\n')
else
    TheEnd = @(VVV) norm(h_I(VVV))<tol;
    fprintf('   Terminator: h_I\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% condiciones iniciales de estimaciOn para todos los algtmos.
V0   = V*eye(N,N);
R0   = V0;
PSI0 = zeros(N,N);
uw0 = [real(V0)*unos; imag(V0)*unos];
rp0 = [R0*unos; PSI0*unos];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P2. ALGORITMO 1: NR-clAsico
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
    %GrAfico de voltajes segUn NR clAsico
    ancho = 200;
    alto  = 200;
    Fontsize = 12;
    
    figure(3)
    subplot(1,2,1)
        plot(0:N,[1;abs(ve_1)/V],'-*')
            box ON
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            title('|v| [p.u.]')
            xlabel('bus No.');
            xlim([0,N]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
    subplot(1,2,2)
        plot(0:N,[0;angle(ve_1)*180/pi],'-*')
            box ON
            grid ON
                hAx=gca;  % avoid repetitive function calls
                set(hAx,'xminorgrid','off','yminorgrid','off')
            title('\angle v [deg.]')
            xlabel('bus No.');
            xlim([0,N]);
            set(gca,'FontName','Times New Roman','FontSize',Fontsize);
            
    set(gcf,'PaperUnits','centimeters',...
            'PaperSize',[2*ancho 1*alto],...
            'PaperPosition',[0 0 2*ancho 1*alto]); %[0 0 ancho alto]
    print('-depsc','-r200','v_NRC') % FunciOn para guardar .eps 
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
    
    fprintf('   No. iteraciones: %i\n',cont)
    fprintf('   Tiempo: %1.1f [ms]\n',elapsedTime*1000)
    fprintf('   ||h_I||: %1.4e\n',norm(h_I(VV)))
    output = sprintf('& %1.1f & %i & %1.2e \n', elapsedTime*1000, cont, norm(h_I(VV)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P3: otros algoritmos y comparar
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
    
    fprintf('   No. iteraciones: %i\n',cont)
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
    
    fprintf('   No. iteraciones: %i\n',cont)
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
    
    fprintf('   No. iteraciones: %i\n',cont)
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
    
    fprintf('   No. iteraciones: %i\n',cont)
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
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%imprimir output
disp(output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
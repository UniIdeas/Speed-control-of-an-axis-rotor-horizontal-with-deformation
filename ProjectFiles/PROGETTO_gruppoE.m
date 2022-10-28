clear all, clc
%% Parameters
mi=0.1; %[kg] Massa concentrata dell'imperferzione
beta=2;
ei=0.5; % distanza massa-asse di rotazione
g=9.8; % accelerazione gravitazionale
Ie=0.75; % momento di inerzia delle eliche senza deformazione
 % ! GARANTITA ROBUSTEZZA PER +- 0.1 !
 
% EQUILIBRI DEL SISTEMA:
teta_e = 0; % accelerazione angolare d'equilibrio
omega_e = 0; %velocità angolare d'equilibrio
u_e = 0; % ingresso di equilibrio
% tau=u(t); coppia applicata al rotore(ingresso del sistema)
  
%RANGE DIAGRAMMA DI BODE:
omega_plot_min = 1e-4;
omega_plot_max = 1e4; 


%(EQUAZIONI DI STATO) System dynamics:
% dot{x_1}=x_2; 
%(mi*ei^2+Ie)*dot{x_2}=-(beta*x_2)-(g*mi*ei*sin(x_1))+teta;
% y=x_2;

%MATRICI SISTEMA LINEARIZZATO:
A=[0,1; -(g*mi*ei)/((mi*ei*ei)+Ie), -beta/((mi*ei*ei)+Ie)];
B=[0;1/((mi*ei*ei)+Ie)];
C=[0,1];
D=0;

%FUNZIONE DI TRASFERIMENTO:
s=tf('s');                                                           
I_s=[s,0;0,s];
GG=2000*s/(1550*s^2+4000*s+981); % GG=C*(inv(I_s-A))*B+D;
GG=zpk(GG); %fattorizzazione
disp("Funzione di Trasferimento G(s):");
display(GG);
grid on, zoom on, hold on;

%FIGURA 1 ( Diagramma di Bode e Margine di fase della G(s) ):
figure(1)
bode(GG) %diagramma di Bode della G(s)
margin(GG)

% Requirements/specifications
W=5;
mu_s=1;
%REGOLATORE STATICO:
Rs = mu_s/(s*s); %da specifica errore a infinito nullo per ingresso (due poli nell'origine)
%w(t)=W*1(t) -> Inserisco un polo nell'origine tramite Rs(s)

% SISTEMA ESTESO:
Ge = Rs*GG;
disp("Sistema esteso");
display(Ge);

%DIAGRAMMA DI BODE DELLA Ge(s):
figure(2)
bode(Ge) 
margin(Ge)

% -> Allowed overshoot S%<=1%=0.01
S_100= 0.01;% sovraelongazione percentuale massima
xi_star= sqrt(log(S_100)^2/(pi^2+log(S_100)^2)); %=0.8261
%margine di fase (rispetta la specifica data: Mf>=45):
Mf=xi_star*100; %=82.6085
display(Mf)

% -> Settling time...Ta,1<=4 [s]:
T_a1Max=4;% tempo di assestamento all' h per cento(con h=1)

%T_a1=4.6/(xi_star*omega_c)<4
%Mf=xi_star*100;
%omega_c>=460/(4*Mf)

omega_cMin=460/(T_a1Max*Mf);

%Punto 4: attenuazione disturbo An volte 
%MAPPATURA SPECIFICHE RUMORE DI MISURA

% Measure Noise requirements 
omega_n=100;
A_n=20*log10(30); %dalla specifica sull'abbattimento del rumore di misura n(t) per omega>100

figure(3);
title("Mappature specifiche e diagramma Ge del sistema esteso");
zoom on, grid on, hold on;
x1 = [omega_n omega_plot_max omega_plot_max omega_n]; 
y1 = [-A_n -A_n 100 100]; 
patch(x1, y1, 'red', 'FaceAlpha', 0.3, 'EdgeAlpha' ,0); %strumento di disegno zona proibita specifiche su n(t)

%Punto 3: omega_cMin zona di attraversamento per L(jw) = 0dB
x21 = [omega_plot_min omega_cMin omega_cMin omega_plot_min];
y21 = [-200 -200 0 0];
patch(x21, y21, 'red', 'FaceAlpha', 0.3, 'EdgeAlpha' ,0); %vincolo su omega_cMin: zona proibita sx

%In mancanza di specifiche imponiamo Wc,max uguale a Wn
x22 = [omega_n omega_plot_max omega_plot_max omega_n];
y22 = [0 0 100 100];
patch(x22, y22, 'red', 'FaceAlpha', 0.5, 'EdgeAlpha' ,0); %vincolo su omega_cMax: zona proibita dx

bodeplot(Ge,{omega_plot_min, omega_plot_max});
margin(Ge);

%Punto 3: MARGINE DI FASE
x3 = [omega_cMin omega_n omega_n omega_cMin];
y3 = [-360 -360 -180+Mf -180+Mf];
patch(x3, y3, 'g', 'FaceAlpha', 0.3, 'EdgeAlpha' ,0);

%Vincolo su margine di fase: 
%il grafico non deve attraversare il range in corrisponde di valori
%precedenti ad L(jw)=0 ma in Wc in quanto Mf è definito in quel punto


%% REGOLATORE DINAMICO
% SCENARIO TIPO B

Rdz = (1+(1/(omega_cMin*1.4)*s)); %aggiungiamo uno zero


figure(4)
mu_d=1.2;
omegac_star=2.5; %Scegliamo omegaC_star nel range [omegac_min , omega_n]
[mag_omega_c_star, arg_omega_c_star, omegac_star] = bode(Ge*Rdz, omegac_star);
amp_Ge_db=20*log10(mag_omega_c_star); 
M_star=10^(-amp_Ge_db/20);


Mf_star=90;% vincolo minimo >=82,60 Mf_star=180+arg(L(jw_c_star));
phi_star=Mf_star-180-arg_omega_c_star;
tau=(M_star- cosd(phi_star))/(omegac_star*sind(phi_star));
alpha_tau=(cosd(phi_star)- 1/M_star)/((omegac_star*sind(phi_star)));

display(M_star);
display(cosd(phi_star)); %vincolo su cos(phi_star)>1/M_star rispettato
display(1/M_star);

%RETE ANTICIPATRICE:
Rds_ant=(1+tau*s)/(1+alpha_tau*s); 

zoom on, grid on, hold on;
x1 = [omega_n omega_plot_max omega_plot_max omega_n]; 
y1 = [-A_n -A_n 100 100];
patch(x1, y1, 'red', 'FaceAlpha', 0.3, 'EdgeAlpha' ,0); %strumento di disegno zona proibita specifiche su n(t)

%Punto 3: omega_cMin zona di attraversamento per L(jw) = 0dB
x21 = [omega_plot_min omega_cMin omega_cMin omega_plot_min];
y21 = [-200 -200 0 0];
patch(x21, y21, 'red', 'FaceAlpha', 0.3, 'EdgeAlpha' ,0); %vincolo su omega_cMin: zona proibita sx

%In mancanza di specifiche imponiamo Wc,max uguale a Wn
x22 = [omega_n omega_plot_max omega_plot_max omega_n];
y22 = [0 0 100 100];
patch(x22, y22, 'red', 'FaceAlpha', 0.5, 'EdgeAlpha' ,0); %vincolo su omega_cMax: zona proibita dx

%FUNZIONE AD ANELLO APERTO:
LS=Rds_ant*Rdz*Ge*mu_d;
bodeplot(LS,{omega_plot_min, omega_plot_max});
margin(LS);


%(Punto 3) vincolo su margine di fase:
x3 = [omega_cMin omega_n omega_n omega_cMin];
y3 = [-360 -360 -180+Mf -180+Mf];
patch(x3, y3, 'g', 'FaceAlpha', 0.3, 'EdgeAlpha' ,0);

%FUNZIONE DI SENSITIVITA' COMPLEMENTARE:
figure(5)
k=1.5;
FF = k*LS/(1+k*LS); %funzione di sensitività complementare
disp("Funzione di sensitività complementare");
display(FF);
display(zpk(FF));
title("BODE PLOT OF F(s)");
bode(FF);

%RISPOSTA AL GRADINO DELLA F(s)
figure(6);
title("Risposta al gradino della F");
%grafico con vincoli su tempo di assestamento e sovraelongazione della
%risposta al gradino della F
hold on;
grid on, zoom on;
x1 = [4 10 10 4];
y1 = [0 0 4.95 4.95];

x2 = [4 10 10 4];
y2 = [5.05 5.05 8 8];

x3 = [0 3.98 3.98 0];
y3 = [5.05 5.05 8 8];

patch(x1, y1, 'y'); %SPECIFICA TEMPO DI ASSESTAMENTO
patch(x2, y2, 'y'); %SPECIFICA TEMPO DI ASSESTAMENTO
patch(x3, y3, 'green'); %SPECIFICA SOVRAELONGAZIONE
step(W*FF, 8);

% Diagramma di bode della S(s)
SS = 1/(1+LS);
figure(7);
title("BODE PLOT OF S(s)");
bode(SS);

%Diagramma di bode della R(s)
RR = mu_d*Rds_ant*Rdz*Rs;
figure(8);
title("BODE PLOT OF R(s)");
bode(RR);
disp("Regolatore");
display(zpk(RR));



%% ANIMAZIONE

[y_step, t_step] = step(W*FF, 4.5);
pause('on');
T=t_step;  
Y=y_step;

fhand=figure(10);
clf;
max_T=max(T);
min_Y=min(Y);
max_Y=max(Y);
Y_ss=Y(end);

for i=1:1:length(T)  

    figure(10);
    clf;
    palag=patch([-0.05 -0.05 0.05 0.05] ,[-ei ei ei -ei],'g'); 
    palab=patch([-ei ei ei -ei] ,[-0.05 -0.05 0.05 0.05],'b'); 
    
    xlim([-1 1])
    ylim([-1 1])

 
    INTEGRANDA = @(t) interp1(T, Y, t);
    integrale = integral(INTEGRANDA,0,T(i));


    dim = [0.45 0.45 0.4 0.4];
    tempi = sprintf('Time = %d s', T(i));
    velocita = sprintf('Omega = %d rad/s', Y(i)+omega_e);
    label = {tempi, velocita};
    annotation('textbox',dim,'String',label,'FitBoxToText','on', 'BackgroundColor','white');


    rotate(palag, [0 0 1], integrale/(2*pi)*360);
    rotate(palab, [0 0 1], integrale/(2*pi)*360);

    pause(0.000001);

if ishandle(fhand)==false; break;
end
end 

%% SIMULINK 
open("SIMULAZIONE_gruppoE"); %apertura file simulink

%Definizione coefficienti del controllore
R = mu_d*Rds_ant*Rdz*Rs;
omega_n=100;
[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};
W=5;
%condizioni iniziali
x0=[0;0];
x_e=[0;0];
u_e=0; %u_e=g*mi*ei*sen(theta_e)
ye=0; %ye = omega_e = 0;
A=[0,1; -981/1550, -80/31];
B=[0;40/31];
C=[0,1];
D=0;
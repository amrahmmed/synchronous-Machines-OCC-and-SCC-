Sn = 2500;         VL = 6600;       f = 50;
polePairs = 1;     Ra = 0.07;       Xl = 0.9;
Nf = 100;          Pcore = 12;      Pfw = 20;             
V = VL/sqrt(3);
Ifl = Sn*1000/(sqrt(3)*VL);
ns = 60*f/polePairs;
w = 2*pi*ns/60;
c = 4;
q = 60/(2* polePairs *3);
alpha = 60/q;
kw = sind(q*alpha/2)/(q*sind(alpha/2));

PF = 0.9; Loading = 0.5;  % 1 Full Load : 0.5 Half Load
phi = -1 *acosd(PF); % 1 = lagging : -1 = leading
Ia = Loading *Ifl;
S = Loading *Sn;  

% OCC
Mdata = [5000 10000 15000 20000 25000 30000];
Edata = [1500 2800 3700 4200 4400 4600];
M = 5000:30000;
E = spline(Mdata,Edata,M);
M = [0 M];
E = [0 E];
k = 1500/5000;        
M0 = interp1(E,M,V);
kmod = V/M0;          
If0 = 2*M0/Nf;       

% SCC at rated armature current
Ersc = Ifl*Xl;
Mrsc = Ersc/k;
Masc = 1.35*c*q*Ifl*kw;
Mfsc = Mrsc + Masc;
Ifsc = 2*Mfsc/Nf;
Isc0 = Ifl*If0/Ifsc;

% Method 1: OCC and SCC 
Xunsat1 = k*M0/Isc0;
Xsat1 = V/Isc0;

% Method 2: short circuit ratio
Ifu = 2*(V/k)/Nf;
SCRu = Ifu/Ifsc;
SCRs = If0/Ifsc;
Zbase = V/Ifl;
Xunsat2 = Zbase/SCRu;
Xsat2 = Zbase/SCRs;

% Exact method
I = Ia*(cosd(phi)-1i*sind(phi));
Er = V + I*(Ra+1i*Xl);
Mr = interp1(E,M,abs(Er));
Ma = 1.35*c*q*Ia*kw;
beta = phi + angle(Er)*180/pi;  
Mf = sqrt(Mr^2 + Ma^2 + 2*Mr*Ma*sind(beta));
exact_IF = 2*Mf/Nf;
exact_Eo = sqrt(3)*interp1(M,E,Mf);
Exact_Voltage_Regulation = 100*(exact_Eo -VL)/VL;

% Approximate method
Ef = abs(V + I*(Ra+1i*Xsat1));
approx_If = 2*Ef/(Nf*kmod);
approximated_E0 = sqrt(3)*Ef; 
approximate_Voltage_Regulation = 100*(approximated_E0 -VL)/VL;

% Output power and prime mover
P = S * PF;                     
Q = S*sind(phi);              
Pcu = 3*Ia^2*Ra/1000;         
Pm = P + Pcu + Pcore + Pfw;   
Torque = 1000*Pm/w;           

disp(['Full Load | 0.9 |leading '])
disp('Xunsat, Xsat . Rows: OCC/SCC, SCR')
disp([Xunsat1 Xsat1; Xunsat2 Xsat2])
disp('If at rated no load voltage, If at rated SC current')
disp([If0 Ifsc])
disp('Rows: exact, approximate. Columns: If(A), E0 line(V), VR(%)')
disp([exact_IF exact_Eo Exact_Voltage_Regulation; approx_If approximated_E0  approximate_Voltage_Regulation])
disp('P(kW),            Q(kvar),         S(kVA),          Pm(kW),        Torque(N m)')
disp([P Q S Pm Torque])

% OCC, air-gap lines and SCC on one plot
If = 2*M/Nf;
Isc = Ifl*(Nf*If/2)/Mfsc;
figure
yyaxis left
plot(If,E,'b',If,k*M,'r--',If,kmod*M,'g-.')
ylabel('Open-circuit phase voltage (V)')
ylim([0 10000])
yyaxis right
plot(If,Isc,'k')
ylabel('Short-circuit current (A)')
ylim([0 800])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
xlabel('Field current (A)')
title('OCC and SCC')
legend('OCC','Air-gap line','New air-gap line','SCC')
grid on
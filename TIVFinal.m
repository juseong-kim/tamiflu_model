dt = 0.01;
time = 0:dt:7;
pV = 210;
beta = 5e-7;
betap = 3e-8;
V0 = 1e+4;
T0 = 7e+7;
gT = 0.8;
deltaV = 5;
deltaI = 2;

ka = 11.04;
ke = 2.64;
EC50 = 30;
emax = 0.98;
omega = 4.63;

V = zeros(1,length(time));
V(1) = V0;
T = V;
T(1) = 7e+7;
I = T;
I(1) = 0;
De = I;
D = I;

%Dadmin = 150; 150 mg tablet 
Dadmin= 75; %75 mg tablet
Td0 = 28/24; %start of treatment day
Td = 12/24; %period of dosage (1 dose per 12 hours)
%Td= 24/24 (one dose per 24 hours)
Td_ind = round((Td0:Td:time(end))./dt+1);
Tadmin = zeros(size(De));
Tadmin(Td_ind) = Dadmin;
%Tadmin(2)= 0; forgetting second dose 


init = [V(1),T(1),I(1),De(1),D(1)]';


options = odeset('RelTol',1e-3,'Abstol',1e-6);




for i = 2:length(time)
    [~,Y] = ode15s(@TIVmodel,[0 dt],init,options,gT,pV,beta,betap,deltaV,deltaI,T0,ka,ke,EC50,emax,omega);
    V(i) = Y(end,1);
    De(i) = Y(end,1);
    T(i) = Y(end,2);
    I(i) = Y(end,3);
    De(i) = Y(end,4) + Tadmin(i);
    D(i) = Y(end,5);
    init = [V(i),T(i),I(i),De(i),D(i)]';
end

clf;

figure();
plot (time, V, 'r');
hold on 
plot (time, T, 'g');
plot (time, I, 'b');
hold off
legend('V ([uV]/cell)', 'T (no. of cells)', 'I (no. of cells)');
xlabel('Time (days)');
title ('TIV model: 75 mg twice/day (forget 1 dose)')



function ynew = TIVmodel(~,y,gT,pV,beta,betap,deltaV,deltaI,T0,ka,ke,EC50,emax,omega)

ynew = zeros(5,1);
ynew(1) = (1-emax*y(5)/(y(5)+ EC50))*pV*y(3)-deltaV*y(1)-beta*y(1)*y(2); %dV/dt
ynew(2) = gT*y(2)*(1-y(2)+y(3))/T0-betap*y(1)*y(2); %dT/dt
ynew(3) = betap*y(1)*y(2)-deltaI*y(3); %dI/dt
ynew(4) = ka*y(4);
ynew(5) = omega*ka*y(4)-ke*y(5); %dD/dt
end 





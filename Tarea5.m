clear 
clc 
close all

syms t s L k

J=0.00048115; %Kg*m^2
Kt=0.22076;   %N*m/A  
Ka=0.22076;   %V/rad*s
Bo=0.0026829;  %N*m*s
Rf=4.08;      %Ohms
Lf=0.011307;  %H

% ---- datos de prueba ----$
% J=1;    %Kg*m^2
% Kt=2;   %N*m/A  
% Ka=2;   %V/rad*s
% Bo=1;   %N*m*s
% Rf=4;   %Ohms
% Lf=2;   %H

num=[Kt];
den=[J*Lf (Bo*Lf+J*Rf) Bo*Rf Ka*Kt];
g=tf(num,den);

Frmx=bandwidth(g);
T=5*Frmx;
To=1/T;

[A B C D]=tf2ss(num,den);

suma=0;
M=3;

%----- discretizacion aproximada -----%
i=A*inv(A);
Ai=i-To*A;
Ad=inv(Ai);
Bd=To*Ad*B;

sys=ss(Ad,Bd,C,D,To);
Gz=tf(sys);

%------- discretizacion exacta -------%
I=s*i;
Ai1=I-A;
Aj=inv(Ai1);
Phi=ilaplace(Aj);
Bi1=Phi*B;
gamma=int(Bi1,0,t);

Ad1=double(subs(Phi,t,To));
Bd1=double(subs(gamma,t,To));

sys1=ss(Ad1,Bd1,C,D,To);
Gz1=tf(sys1);

%------ discretizacion truncada ------%
for n=1:1:M 
    s = ((A^n*To^n)/(n)) + suma;
    suma=s;
    Phi2=i+s;
    Bi2=Phi2*B;
    gamma2=int(Bi2,0,t);
end

Bd2=double(subs(gamma2,t,To));

sys2=ss(Phi2,Bd2,C,D,To);
Gz2=tf(sys2);

figure(1)
step(g)
hold on 
step(Gz)
title('Representacion del sistema continua vs discretizado aproximado');
legend('continua','discreta');
grid on

figure(2)
step(g)
hold on 
step(Gz1)
title('Representacion del sistema continua vs discretizado Exacta');
legend('continua','discreta');
grid on

figure(3)
step(g)
hold on 
step(Gz2)
title('Representacion del sistema continua vs discretizado Truncada');
legend('continua','discreta');
grid on


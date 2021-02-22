%% Laboratorio 3
%% Datos
% Hector Alejandro Klée González - 17118
% Juan Diego Castillo Amaya - 17074
%% Función de Transferencia Simbólica
% Definición de Variables
syms Vout Vin s r1_sym r2_sym r3_sym c1_sym c2_sym c3_sym
% Simplificaciones
% Nodo C
Vc = Vout;
Vb = Vc * ((r3_sym * c2_sym * s) + 1); 
% Nodo B
Va = Vb + (((((Vb - Vout) * c3_sym * s)) - ((Vc - Vb) / r3_sym)) * r2_sym);
% Nodo A
Vin = Va + (((Va * c1_sym * s) + ((Va - Vb) / r2_sym)) * r1_sym);
Vin = simplify(Vin);
% Función de Transferencia Simbólica
G_sym = Vout / Vin;
G_sym = simplify(G_sym);
%% Primera Parte
% Definición de Variables
b0_sym = 1;
a0_sym = 1;
a1_sym = (c1_sym * r1_sym) + (c2_sym * r1_sym) + (c2_sym * r2_sym) + ...
    (c2_sym * r3_sym);
a2_sym = (c1_sym * c2_sym * r1_sym * r2_sym)...
    +(c1_sym * c2_sym * r1_sym * r3_sym)...
    +(c2_sym * c3_sym * r1_sym * r3_sym)...
    +(c2_sym * c3_sym * r2_sym * r3_sym);
a3_sym = (c1_sym * c2_sym * c3_sym * r1_sym * r2_sym * r3_sym);

% Realización Controlable
p1_RC = -[a0_sym a1_sym a2_sym];
p2_RC = [zeros(length(p1_RC) - 1, 1) eye(length(p1_RC) - 1)];
A1_sym = [p2_RC; p1_RC];
B1_sym = [0; 0; 1];
C1_sym = [b0_sym 0 0];

% Realización Observable
p1_RO = -[a2_sym; a1_sym; a0_sym];
p2_RO = [eye(length(p1_RO) - 1); zeros(1, length(p1_RO) - 1)];
A2_sym = [p1_RO p2_RO];
B2_sym = [0; 0; b0_sym];
C2_sym = [1 0 0];

% Ecuación de Diferencias

%% Segunda Parte
%% Función de Transferencia Numérica
% Definición de Variables
syms Vout Vin s 
% Definición de Variables
r1 = 1e3;
r2 = 10e3;
r3 = r2;
c1 = 1e-6;
c2 = 0.1e-6;
c3 = 10e-6;
% Simplificaciones
% Nodo C
Vc = Vout;
Vb = Vc * ((r3 * c2 * s) + 1); 
% Nodo B
Va = Vb + (((((Vb - Vout) * c3 * s)) - ((Vc - Vb) / r3)) * r2);
% Nodo A
Vin = Va + (((Va * c1 * s) + ((Va - Vb) / r2)) * r1);
Vin = simplify(Vin);
% Función de Transferencia Simbólica
G0_sym = Vout / Vin;
G0_sym = simplify(G0_sym);

% Función de Transferencia tf
% Obtenención de numerador y denominador
[N, D] = numden(simplifyFraction(G0_sym));
% Obtenención coeficientes del numerador y denominador
b = fliplr(double(coeffs(N)));
a = fliplr(double(coeffs(D)));
% Expresión tipo tf
G0 = tf(b, a);

% Inciso 1
% Definición de Variables
b0 = G0.Numerator{1}(4);
a0 = G0.Denominator{1}(4);
a1 = G0.Denominator{1}(3);
a2 = G0.Denominator{1}(2);
a3 = G0.Denominator{1}(1);

% Realización Controlable
p1_RC_num = -[a0 a1 a2];
p2_RC_num = [zeros(length(p1_RC_num) - 1, 1) eye(length(p1_RC_num) - 1)];
A1 = [p2_RC_num; p1_RC_num];
B1 = [0; 0; 1];
C1 = [G0.Numerator{1}(4) 0 0];

% Realización Observable
p1_RO_num = -[a2; a1; a0];
p2_RO_num = [eye(length(p1_RO_num) - 1); zeros(1, length(p1_RO_num) - 1)];
A2 = [p1_RO_num p2_RO_num];
B2 = [0; 0; G0.Numerator{1}(4)];
C2 = [1 0 0];

A3=[-((r2+r1)/(r1*r2*c1)),1/(r2*c1),-1/(r2*c1);...
    0,0,-1/(r3*c2);...
    -1/(r2*c3),1/(r2*c3),-(r3+r2)/(r2*r3*c3)];
B3=[1/(r1*c1);0;0];
C3=[-0 1 0];
% Inciso 3
[A4,B4,C4,D4] = tf2ss(b,a);

% Inciso 4
% save linSys.mat linsys1
load linSys.mat linsys1
A5 = linsys1.A;
B5 = linsys1.B;
C5 = linsys1.C;
D5 = linsys1.D;

% Inciso 5
s = tf('s');
G1 = C1 * (inv((s * eye(3)) - A1)) * B1;
G2 = C2 * (inv((s * eye(3)) - A2)) * B2;
G3 = C3 * (inv((s * eye(3)) - A3)) * B3;
G4 = C4 * (inv((s * eye(3)) - A4)) * B4;
G5 = C5 * (inv((s * eye(3)) - A5)) * B5;
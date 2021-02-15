%% Laboratorio 3
%% Datos
% Hector Alejandro Kl�e Gonz�lez - 17118
% Juan Diego Castillo Amaya - 17074
%% Funci�n de Transferencia Simb�lica
% Definici�n de Variables
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
% Funci�n de Transferencia Simb�lica
G_sym = Vout / Vin;
G_sym = simplify(G_sym);
%% Primera Parte
% Definici�n de Variables
b0_sym = 1;
a0_sym = 1;
a1_sym = (c1_sym * r1_sym) + (c2_sym * r1_sym) + (c2_sym * r2_sym) + ...
    (c2_sym * r3_sym);
a2_sym = (c1_sym * c2_sym * r1_sym * r2_sym)...
    +(c1_sym * c2_sym * r1_sym * r3_sym)...
    +(c2_sym * c3_sym * r1_sym * r3_sym)...
    +(c2_sym * c3_sym * r2_sym * r3_sym);
a3_sym = (c1_sym * c2_sym * c3_sym * r1_sym * r2_sym * r3_sym);

% Realizaci�n Controlable
p1_RC = -[a0_sym a1_sym a2_sym a3_sym];
p2_RC = [zeros(length(p1_RC) - 1, 1) eye(length(p1_RC) - 1)];
A1_sym = [p2_RC; p1_RC];
B1_sym = 1;
C1_sym = b0_sym;

% Realizaci�n Observable
p1_RO = -[a3_sym; a2_sym; a1_sym; a0_sym];
p2_RO = [eye(length(p1_RO) - 1); zeros(1, length(p1_RO) - 1)];
A2_sym = [p1_RO p2_RO];
B2_sym = b0_sym;
C2_sym = 1;

% Ecuaci�n de Diferencias

%% Segunda Parte
%% Funci�n de Transferencia Num�rica
% Definici�n de Variables
syms Vout Vin s 
% Definici�n de Variables
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
% Funci�n de Transferencia Simb�lica
G0_sym = Vout / Vin;
G0_sym = simplify(G0_sym);

% Funci�n de Transferencia tf
% Obtenenci�n de numerador y denominador
[N, D] = numden(simplifyFraction(G0_sym));
% Obtenenci�n coeficientes del numerador y denominador
b = fliplr(double(coeffs(N)));
a = fliplr(double(coeffs(D)));
% Expresi�n tipo tf
G0 = tf(b, a);

% Inciso 1
% Definici�n de Variables
b0 = G0.Numerator{1}(4);
a0 = G0.Denominator{1}(1);
a1 = G0.Denominator{1}(2);
a2 = G0.Denominator{1}(3);
a3 = G0.Denominator{1}(4);

% Realizaci�n Controlable
p1_RC_num = fliplr(-G0.Denominator{1});
p2_RC_num = [zeros(length(p1_RC_num) - 1, 1) eye(length(p1_RC_num) - 1)];
A1 = [p2_RC_num; p1_RC_num];
B1 = 1;
C1 = G0.Numerator{1}(4);

% Realizaci�n Observable
p1_RO_num = (-G0.Denominator{1})';
p2_RO_num = [eye(length(p1_RO_num) - 1); zeros(1, length(p1_RO_num) - 1)];
A2 = [p1_RO_num p2_RO_num];
B2 = G0.Numerator{1}(4);
C2 = 1;

% Inciso 3
[A4,B4,C4,D4] = tf2ss(b,a);
%% Laboratorio 3
%% Datos
% Hector Alejandro Kl�e Gonz�lez - 17118
% Juan Diego Castillo Amaya - 17074
%% Funci�n de Transferencia Simb�lica
% Definici�n de Variables
syms Vout Vin s r1 r2 r3 c1 c2 c3
% Definici�n de Variables
%{
r1 = 1e3;
r2 = 10e3;
r3 = r2;
c1 = 1e-6;
c2 = 0.1e-6;
c3 = 10e-6;
%}
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
G_sym = Vout / Vin;
G_sym = simplify(G_sym);
%% Funci�n de Transferencia tf
% Obtenenci�n de numerador y denominador
[N, D] = numden(simplifyFraction(G_sym));
% Obtenenci�n coeficientes del numerador y denominador
b = fliplr(double(coeffs(N)));
a = fliplr(double(coeffs(D)));
% Expresi�n tipo tf
G0 = tf(b, a);
%% Primera Parte
% Definici�n de Variables
b0 = 1;
a0 = 1;
a1 = (c1 * r1) + (c2 * r1) + (c2 * r2) + (c2 * r3);
a2 = (c1 * c2 * r1 * r2) + (c1 * c2 * r1 * r3) + ...
    (c2 * c3 * r1 * r3) + (c2 * c3 * r2 * r3);
a3 = (c1 * c2 * c3 * r1 * r2 * r3);
% Realizaci�n Controlable
p1_RC = -[a0 a1 a2 a3];
p2_RC = [zeros(length(p1_RC) - 1, 1) eye(length(p1_RC) - 1)];
A1_sym = [p2_RC ;p1_RC];
B1_sym = 1;
C1_sym = b0;

% Realizaci�n Observable
p1_RO = -[a3; a2; a1; a0];
p2_RO = [eye(length(p1_RO) - 1); zeros(1, length(p1_RO) - 1)];
A2_sym = [p1_RO p2_RO];
B2_sym = b0;
C2_sym = 1;
% Ecuaci�n de Diferencias
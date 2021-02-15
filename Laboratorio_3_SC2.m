%% Laboratorio 3
%% Datos
% Hector Alejandro Kl�e Gonz�lez - 17118
% Juan Diego Castillo Amaya - 17074
%% Funci�n de Transferencia Simb�lica
% Definici�n de Variables
syms Vout Vin s %r1 r2 r3 c1 c2 c3
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
%Vb = simplify(Vb);
% Nodo B
Va = Vb + (((((Vb - Vout) * c3 * s)) - ((Vc - Vb) / r3)) * r2);
%Va = simplify(Va);
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
% Normalizamos
%b = b / a(1);
%a = a / a(1);
G = tf(b, a);
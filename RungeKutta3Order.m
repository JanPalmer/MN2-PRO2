function [Y, X] = RungeKutta3Order(a, b, n, y0, f, p, q, r, dy0)
% Autor: Jan Palmer 313429
%
% Funkcja RungeKutta3Order zwraca dwa wektory poziome: wektor 'X',
% zawierajacy n+1 punktow rownomiernie rozlozonych na przedziale [a, b],
% oraz wektor 'Y', zawierajacy przyblizenia funkcji opisanej rozpatrywanym
% rownaniem rozniczkowym.
% Kazdemu z punktow w 'X' odpowiada wartosc w wektorze 'Y', przyblizona
% przy pomocy metody Rungego Kutty rzedu 3-go, na podstawie podanych
% parametrow liniowego rownania rozniczkowego 1-go lub 2-go rzedu.
% Dla rownan rozniczkowych 2-go rzedu, wektor 'Y' bedzie wektorem
% dwurzedowym, gdzie 1. rzad to przyblizenia y(x), a 2. to przyblizenia
% y'(x).
%
% Przyjmujemy, ze rownanie moze wystepowac w nastepujacych postaciach:
% q(x)y' + p(x)y = f(x), jako rownanie rozniczkowe 1-go rzedu,
% lub
% r(x)y'' + q(x)y' + p(x)y = f(x), jako rownanie rozniczkowe 2-go rzedu.
%
% Ze wzgledu na wzor F(x, Y_k), dla rownan rozniczkowych 1-go rzedu
% funkcja q(x) nie moze osiagac wartosci 0 na przedziale [a, b]. To samo
% dotyczy funkcji r(x) dla rownan rozniczkowych 2-go rzedu.
%
% WEJSCIE:
%       a - poczatek przedzialu, na ktorym przyblizac bedziemy badane
%           rownanie rozniczkowe.
%       b - koniec rozpatrywanego przedzialu.
%       n - liczba podprzedzialow, na ktore podzielony zostanie
%           przedzial [a, b].
%      y0 - warunek poczatkowy y(a).
%       f - uchwyt do funkcji znajdujacej sie po prawej stronie rownania
%           rozniczkowego
%       p - uchwyt do funkcji znajdujacej sie przy y(x).
%       q - uchwyt do funkcji znajdujacej sie przy y'(x). Funkcja nie moze
%           osiagac wartosci 0 dla rownan rozniczkowych 1-go rzedu.
%       r - argument opcjonalny, wykorzystywany dla rownan rozniczkowych
%           2-go rzedu. Uchwyt do funkcji znajdujacej sie przy y''(x).
%           Funkcja nie moze osiagac wartosci 0.
%     dy0 - argument opcjonalny, wykorzystywany dla rownan rozniczkowych
%           2-go rzedu. Warunek poczatkowy y'(a).
% WYJSCIE:
%       Y - wektor przyblizen rozwiazan rownania rozniczkowego. Dla rownan
%           2-go rzedu, Y bedzie 2-rzedowym wektorem, gdzie 1 rzad.
%           to przyblizenia y(x), a 2 rzad to przyblizenia y'(x).
%       X - wektor wytyczonych punktow na osi x, na ktorych obliczone
%           zostalo przyblizenie rozwiazan rownania rozniczkowego

% Sprawdzenie rzedu rownania rozniczkowego npdst podanych argumentow
switch nargin
    case 9
        Is2Order = 1;
    case 7
        Is2Order = 0;
    otherwise
        return;
end % switch

% Stale
h = (b - a)/n;
X = linspace(a, b, n+1);
gamma = [1/2, 0; -1, 2];
c = [1/6, 2/3, 1/6];

% Wyznaczenie rownania funkcji F(x, Y_k) w zaleznosci od rzedu
% rozpatrywanego rownania rozniczkowego
if Is2Order == 0
    Y = zeros(1, n+1);
    Y(1) = y0;
    F = @(x, Y_k) (f(x) - p(x)*Y_k)./q(x);
else
    Y = zeros(2, n+1);
    Y(1, 1) = y0;
    Y(2, 1) = dy0;
    F = @(x, Y_k) [Y_k(2); (f(x) - p(x)*Y_k(1) - q(x)*Y_k(2))./r(x)];
end % if

% Przyblizanie Y przy pomocy metody Rungego-Kutty rzedu 3-go
for i = 1:n
    K0 = F(X(i), Y(:, i));
    K1 = F(X(i) + h*gamma(1, 1), Y(:, i) + h*gamma(1, 1)*K0);
    K2 = F(X(i) + h*(gamma(2, 1) + gamma(2, 2)), Y(:, i) ...
        + h*(gamma(2, 1)*K0 + gamma(2, 2)*K1));
    Y(:, i+1) = Y(:, i) + h*(c(1)*K0 + c(2)*K1 + c(3)*K2);
end % for

end % function
function [Y, X] = AdamsBashforthMain(a, b, n, y0, f, p, q, r, dy0)
% Autor: Jan Palmer 313429
%
% Funkcja AdamsBashforthMain zwraca dwa wektory poziome: wektor 'X',
% zawierajacy n+1 punktow rownomiernie rozlozonych na przedziale [a, b],
% oraz wektor 'Y', zawierajacy przyblizenia funkcji opisanej rozpatrywanym
% rownaniem rozniczkowym.
% Kazdemu z punktow w 'X' odpowiada wartosc w wektorze 'Y', przyblizona
% przy pomocy metody Adamsa-Bashfortha rzedu 3-go, na podstawie podanych
% parametrow liniowego rownania rozniczkowego 1-go lub 2-go rzedu.
% Dla rownan rozniczkowych 2-go rzedu, wektor 'Y' bedzie wektorem
% dwurzedowym, gdzie 1. rzad to przyblizenia y(x), a 2. to przyblizenia
% y'(x).
% Jako ze metoda Adamsa-Bashfortha jest metoda 3-krokowa, pierwsze 2 kroki
% po punkcie wyznaczonym przez warunki poczatkowe zostaly wyznaczone metoda
% Rungego-Kutty 3-go rzedu.
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
%           przedzial [a, b]. Nie moze byc mniejsza niz 2.
%      y0 - warunek poczatkowy y(a).
%       f - uchwyt do funkcji znajdujacej sie po prawej stronie rownania
%           rozniczkowego
%       p - uchwyt do funkcji znajdujacej sie przy y(x).
%       q - uchwyt do funkcji znajdujacej sie przy y'(x). Funkcja nie moze
%           osiagac wartosci 0 dla rownan rozniczkowych 1 rzedu.
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

% Wyznaczenie rownania funkcji F(x, Y_k) w zaleznosci od rzedu
% rozpatrywanego rownania rozniczkowego oraz zainicjowanie pierwszych
% 3 wartosci Y przy pomocy metody Rungego-Kutty rzedu 3-go
if Is2Order == 0
    Y = zeros(1, n+1);
    F = @(x, Y_k) (f(x) - p(x)*Y_k)./q(x);
    Y_RK = RungeKutta3Order(X(1), X(3), 2, y0, f, p, q);
    Y(1:3) = Y_RK(1:3);
else
    Y = zeros(2, n+1);
    F = @(x, Y_k) [Y_k(2); (f(x) - p(x)*Y_k(1) - q(x)*Y_k(2))./r(x)];
    Y_RK = RungeKutta3Order(X(1), X(3), 2, y0, f, p, q, r, dy0);
    Y(:, 1:3) = Y_RK(:, 1:3);
end % if

% Przyblizanie Y przy pomocy metody Adamsa-Bashfortha rzedu 3-go
for i = 3:n
    Y(:, i+1) = Y(:, i) + (h/12)*(23*F(X(i), Y(:, i)) ...
        - 16*F(X(i-1), Y(:, i-1)) + 5*F(X(i-2), Y(:, i-2)));
end % for

end % function
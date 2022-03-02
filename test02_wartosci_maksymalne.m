function test02_wartosci_maksymalne()
% Autor: Jan Palmer 313429
%
% Funkcja porownuje wartosci maksimum i minimum na danym przedziale  
% otrzymanych przy pomocy przyblizenia metoda Rungego-Kutty rzedu 3 
% oraz Adamsa-Bashfortha rzedu 3 z wartoscia uzyskana przy pomocy wzoru 
% funkcji, ktora odpowiada funkcji 'y(x)' opisanej przez badane rownanie 
% rozniczkowe, lub z wartoscia uzyskana przy pomocy procedury 'ode45', 
% z opcjami 'RelTol',2.22045e-14,'AbsTol',1e-16.
% Badana jest takze dokladnosc przyblizenia zaleznie od ilosci 
% podprzedzialow, ktorymi podzielony zostanie przedzial [a, b].
% Jako 'AB' bedzie oznaczana metoda Adamsa-Bashfortha 3-rzedu, a jako 
% 'RK' - metoda Rungego-Kutty 3-rzedu.
% Dla niektorych funkcji maksimum lub minimum znajduje sie w punkcie
% okreslonym przez warunki poczatkowe, dlatego tez dla niektorych 
% tych wartosci blad wynosi 0 niezaleznie od ilosci podprzedzialow.

clearvars
format longe

disp(['Test dokladnosci przyblizen wartosci maksimum i minimum na ' ...
    'rozpatrywanych przedzialach.']);

Err = zeros(4, 3);
disp(' ');
disp("Test rownan rozniczkowych rzedu 1-go");

disp("   1. y'' + 2y' + 2y = 0, y(0) = 0, y'(0) = 5, a = 0, b = 6");
n = [100, 10000];
f = @(x) 0;
p = @(x) 2;
q = @(x) 2;
r = @(x) 1;
y0 = 0;
dy0 = 5;
a = 0;
b = 6;

for i = 1:2
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q, r, dy0);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q, r, dy0);
    Y = 5.*exp(-X).*sin(X);
    Y_max = max(Y);
    Y_min = min(Y);
    Err(1, i) = abs(max(Y_RK(1, :)) - Y_max);
    Err(2, i) = abs(max(Y_AB(1, :)) - Y_max);
    Err(3, i) = abs(min(Y_RK(1, :)) - Y_min);
    Err(4, i) = abs(min(Y_AB(1, :)) - Y_min);
end % for

disp("Blad maksimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2))]);
disp("Blad minimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(3, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(3, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(4, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(4, 2))]);

disp("   2. (x/2)y' + (-1)y = -10, a = -10, b = -2, y(-10) = 4");
f = @(x) -10;
p = @(x) -1;
q = @(x) x/2;
y0 = 4;
a = -10;
b = -2;

for i = 1:2
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q);
    Y = 10 - (3.*X.*X)./50;
    Y_max = max(Y);
    Y_min = min(Y);
    Err(1, i) = abs(max(Y_RK(1, :)) - Y_max);
    Err(2, i) = abs(max(Y_AB(1, :)) - Y_max);
    Err(3, i) = abs(min(Y_RK(1, :)) - Y_min);
    Err(4, i) = abs(min(Y_AB(1, :)) - Y_min);
end % for

disp("Blad maksimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2))]);
disp("Blad minimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(3, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(3, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(4, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(4, 2))]);

disp("   3. (x/3)y' + (-1)y = -10, a = 10, b = 15, y(10) = 4");
f = @(x) -10;
p = @(x) -1;
q = @(x) x/3;
y0 = 4;
a = 10;
b = 15;

for i = 1:2
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q);
    Y = 10 - (3.*X.*X.*X)./500;
    Y_max = max(Y);
    Y_min = min(Y);
    Err(1, i) = abs(max(Y_RK(1, :)) - Y_max);
    Err(2, i) = abs(max(Y_AB(1, :)) - Y_max);
    Err(3, i) = abs(min(Y_RK(1, :)) - Y_min);
    Err(4, i) = abs(min(Y_AB(1, :)) - Y_min);
end % for

disp("Blad maksimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2))]);
disp("Blad minimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(3, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(3, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(4, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(4, 2))]);

disp("   4. (x/3)y' + y = -10, a = -10, b = -8, y(10) = -24");
f = @(x) -10;
p = @(x) 1;
q = @(x) x/3;
y0 = -24;
a = -10;
b = -8;

for i = 1:2
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q);
    Y = 14000./(X.*X.*X) - 10;
    Y_max = max(Y);
    Y_min = min(Y);
    Err(1, i) = abs(max(Y_RK(1, :)) - Y_max);
    Err(2, i) = abs(max(Y_AB(1, :)) - Y_max);
    Err(3, i) = abs(min(Y_RK(1, :)) - Y_min);
    Err(4, i) = abs(min(Y_AB(1, :)) - Y_min);
end % for

disp("Blad maksimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2))]);
disp("Blad minimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(3, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(3, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(4, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(4, 2))]);

disp(' ');
disp("Test rownan rozniczkowych rzedu 2-go");

disp("   1. xy'' + 0y' + 0y = 1, a = 1, b = 13, y(1) = 0, y'(1) = 0");
f = @(x) 1;
p = @(x) 0;
q = @(x) 0;
r = @(x) x;
y0 = 0;
dy0 = 0;
a = 1;
b = 13;

for i = 1:2
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q, r, dy0);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q, r, dy0);
    Y = -X + X.*log(X) + 1;
    Y_max = max(Y);
    Y_min = min(Y);
    Err(1, i) = abs(max(Y_RK(1, :)) - Y_max);
    Err(2, i) = abs(max(Y_AB(1, :)) - Y_max);
    Err(3, i) = abs(min(Y_RK(1, :)) - Y_min);
    Err(4, i) = abs(min(Y_AB(1, :)) - Y_min);
end % for

disp("Blad maksimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2))]);
disp("Blad minimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(3, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(3, 2))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(4, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(4, 2))]);

disp("   2. (-x)y'' + 1/2y' + y = 1, a = 1, b = 100, y(1) = 0, y'(1) = 0");
n = [100, 100000, 500000];
f = @(x) 1;
p = @(x) 1;
q = @(x) 1/2;
r = @(x) -x;
y0 = 0;
dy0 = 0;
a = 1;
b = 100;

dydt = @(t, y) [y(2); (-y(2)/2 + 1 - y(1))./(-t)];
opts = odeset('RelTol',2.22045e-14,'AbsTol',1e-16);
sol = ode45(dydt, [a b+0.0001], [y0; dy0], opts);

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q, r, dy0);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q, r, dy0);
    Y = deval(sol, X);
    Y_max = max(Y(1, :));
    Y_min = min(Y(1, :));
    Err(1, i) = abs(max(Y_RK(1, :)) - Y_max);
    Err(2, i) = abs(max(Y_AB(1, :)) - Y_max);
    Err(3, i) = abs(min(Y_RK(1, :)) - Y_min);
    Err(4, i) = abs(min(Y_AB(1, :)) - Y_min);
end % for

disp("Blad maksimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)) ',' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)) ',' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3))]);
disp("Blad minimum");
disp(['RK - n = ' num2str(n(1)) ', Blad = ' num2str(Err(3, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(3, 2)) ',' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(3, 3))]);
disp(['AB - n = ' num2str(n(1)) ', Blad = ' num2str(Err(4, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(4, 2)) ',' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(4, 3))]);

end % function
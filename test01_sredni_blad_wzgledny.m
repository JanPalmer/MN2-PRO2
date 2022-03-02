function test01_sredni_blad_wzgledny()
% Autor: Jan Palmer 313429
% 
% Funkcja wyznacza sredni blad bezwzgledny dla roznych rownan rozniczkowych 
% 1-go i 2-go rzedu, porownujac wyniki otrzymane przy pomocy metody 
% Rungego-Kutty 3-go rzedu oraz metody Adamsa-Bashfortha 3-go rzedu
% z rzeczywistymi wartosciami funkcji dla rownan 1-go rzedu, oraz z
% przyblizeniami procedura ode45 dla rownan 2-go rzedu (oprocz ostatniego).

Err = zeros(2, 3);
disp(' ');
disp("Test rownan rozniczkowych rzedu 1-go");

disp("   1. 1y' + 1y = 1, a = 0, b = 6, y(0) = 0");
n = [60, 6000, 60000];
f = @(x) 1;
a = 0;
b = 6;
y0 = 0;

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, f, f);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, f, f);
    Y = -exp(-X) + 1;
    Err(1, i) = sum(abs(Y_AB - Y)/Y)/length(X);
    Err(2, i) = sum(abs(Y_RK - Y)/Y)/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3)*100) '%']);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3)*100) '%']);

disp("   2. -0.5y' + 2y = 0, a = 0, b = 5, y(0) = -10");
n = [100, 10000, 100000];
f = @(x) 0;
p = @(x) -1/2;
q = @(x) 2;
y0 = -10;
a = 0;
b = 5;

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q);
    Y = -10.*exp(4.*X);
    Err(1, i) = sum(abs(Y_AB - Y)/Y)/length(X);
    Err(2, i) = sum(abs(Y_RK - Y)/Y)/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3)*100) '%']);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3)*100) '%']);

disp("   3. sin(x)y' + y = -pi, a = -2, b = -0.5, y(-2) = pi/2");
n = [100, 10000, 100000];
f = @(x) pi/2;
p = @(x) 1;
q = @(x) sin(x);
y0 = pi/2;
a = -2;
b = -0.5;

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q);
    Y = -0.5.*pi*(3.*tan(1).*cot(X./2) + 2);
    Err(1, i) = sum(abs(Y_AB - Y)/Y)/length(X);
    Err(2, i) = sum(abs(Y_RK - Y)/Y)/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3)*100) '%']);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3)*100) '%']);

disp("   4. (x/2)y' + (-1)y = -10, a = -10, b = -1, y(-10) = 4");
n = [100, 10000, 100000];
f = @(x) -10;
p = @(x) -1;
q = @(x) x/2;
y0 = 4;
a = -10;
b = -1;

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q);
    Y = 10 - (3.*X.*X)./50;
    Err(1, i) = sum(abs(Y_AB - Y)/Y)/length(X);
    Err(2, i) = sum(abs(Y_RK - Y)/Y)/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3)*100) '%']);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3)*100) '%']);

disp(' ');
disp("Test rownan rozniczkowych rzedu 2-go");

disp("   1. 1/2y'' + (-1/2)y' + 2y= 0, y(-15) = 0, y'(-15) = 0");
n = [10, 1000, 10000];
f = @(x) 1/9;
p = @(x) 2;
q = @(x) -1/2;
r = @(x) 1/2;
y0 = 0;
dy0 = 0;
a = -15;
b = 15;

dydt = @(t, y) [y(2); ((1/2)*y(2) + (-2)*y(1) + 1/9)*2];
opts = odeset('RelTol',2.22045e-14,'AbsTol',1e-16);
sol = ode45(dydt, [a b+0.0001], [y0; dy0], opts);

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q, r, dy0);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q, r, dy0);
    Y = deval(sol, X);
    Err(1, i) = sum(abs(Y_AB(1, :) - Y(1, :))/Y(1, :))/length(X);
    Err(2, i) = sum(abs(Y_RK(1, :) - Y(1, :))/Y(1, :))/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)) ',' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3))]);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)) ','  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)) ',' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3))]);

disp("   2. 8y'' + 4y' + (-2)y = 1/3, a = 100, b = 200, y(100) = 0, y'(100) = 2");
n = [10, 1000, 10000];
f = @(x) 1/3;
p = @(x) -2;
q = @(x) 4;
r = @(x) 8;
y0 = 0;
dy0 = 2;
a = 100;
b = 200;

dydt = @(t, y) [y(2); (-4*y(2) + 2*y(1) + 1/3)/8];
opts = odeset('RelTol',2.22045e-14,'AbsTol',1e-16);
sol = ode45(dydt, [a b+0.0001], [y0; dy0], opts);

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q, r, dy0);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q, r, dy0);
    Y = deval(sol, X);
    Err(1, i) = sum(abs(Y_AB(1, :) - Y(1, :))/Y(1, :))/length(X);
    Err(2, i) = sum(abs(Y_RK(1, :) - Y(1, :))/Y(1, :))/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3)*100) '%']);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3)*100) '%']);

disp("   3. xy'' + (-x)y' + 0y = 1, a = 1, b = 100, y(1) = 0, y'(1) = 0");
n = [10, 1000, 10000];
f = @(x) 1;
p = @(x) 0;
q = @(x) -x;
r = @(x) x;
y0 = 0;
dy0 = 0;
a = 1;
b = 100;

dydt = @(t, y) [y(2); (t*y(2) + 1)/t];
opts = odeset('RelTol',2.22045e-14,'AbsTol',1e-16);
sol = ode45(dydt, [a b+0.0001], [y0; dy0], opts);

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q, r, dy0);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q, r, dy0);
    Y = deval(sol, X);
    Err(1, i) = sum(abs(Y_AB(1, :) - Y(1, :))/Y(1, :))/length(X);
    Err(2, i) = sum(abs(Y_RK(1, :) - Y(1, :))/Y(1, :))/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3)*100) '%']);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3)*100) '%']);

disp("   4. xy'' + 0y' + 0y = 1, a = 1, b = 13, y(1) = 0, y'(1) = 0");
n = [10, 1000, 10000];
f = @(x) 1;
p = @(x) 0;
q = @(x) 0;
r = @(x) x;
y0 = 0;
dy0 = 0;
a = 1;
b = 13;

for i = 1:3
    [Y_AB, X] = AdamsBashforthMain(a, b, n(i), y0, f, p, q, r, dy0);
    Y_RK = RungeKutta3Order(a, b, n(i), y0, f, p, q, r, dy0);
    Y = -X + X.*log(X) + 1;
    Err(1, i) = sum(abs(Y_AB(1, :) - Y(1, :))/Y(1, :))/length(X);
    Err(2, i) = sum(abs(Y_RK(1, :) - Y(1, :))/Y(1, :))/length(X);
end % for

disp("Blad dla metody Rungego-Kutty rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(2, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(2, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(2, 3)*100) '%']);
disp("Blad dla metody Adamsa-Bashfortha rzedu 3-go");
disp(['n = ' num2str(n(1)) ', Blad = ' num2str(Err(1, 1)*100) '%,'  ...
    ' n = ' num2str(n(2)) ', Blad = ' num2str(Err(1, 2)*100) '%,' ...
    ' n = ' num2str(n(3)) ', Blad = ' num2str(Err(1, 3)*100) '%']);

end % function
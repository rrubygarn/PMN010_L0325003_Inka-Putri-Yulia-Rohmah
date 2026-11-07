%L0325003_Inka Putri
%Persamaan dan parameter
f  = @(x, y) (4*x.^3 - 2*x + 1) ./ (2 + y);
x0 = 0;  y0 = 1;  xn = 2;  n = 15;
h  = (xn - x0) / n;

%Solusi analitik (eksak)
% (2+y) dy = (4x^3 - 2x + 1) dx
% 2y + y^2/2 = x^4 - x^2 + x + C, dengan y(0)=1 -> C = 2.5
% y = -2 + sqrt(2x^4 - 2x^2 + 2x + 9)
yex = @(x) -2 + sqrt(2*x.^4 - 2*x.^2 + 2*x + 9);

x = x0 + (0:n) * h;          % titik-titik x
Ye = yex(x);                 % nilai eksak di setiap titik

%Metode Euler
Yeu = zeros(1, n+1);  Yeu(1) = y0;
for i = 1:n
  Yeu(i+1) = Yeu(i) + h * f(x(i), Yeu(i));
end

%Metode Heun
Yhe = zeros(1, n+1);  Yhe(1) = y0;
for i = 1:n
  k1 = f(x(i), Yhe(i));
  yp = Yhe(i) + h * k1;                 % prediktor (Euler)
  k2 = f(x(i+1), yp);
  Yhe(i+1) = Yhe(i) + (h/2) * (k1 + k2); % korektor
end

%Metode Runge-Kutta Orde-4
Yrk = zeros(1, n+1);  Yrk(1) = y0;
for i = 1:n
  k1 = f(x(i),       Yrk(i));
  k2 = f(x(i) + h/2, Yrk(i) + h/2 * k1);
  k3 = f(x(i) + h/2, Yrk(i) + h/2 * k2);
  k4 = f(x(i) + h,   Yrk(i) + h * k3);
  Yrk(i+1) = Yrk(i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

%Galat pada x = 2
y_eksak = Ye(end);
hasil   = [Yeu(end), Yhe(end), Yrk(end)];
nama    = {'Euler', 'Heun', 'Runge-Kutta 4'};

galat_abs = abs(y_eksak - hasil);            % galat abstrak (mutlak)
galat_rel = galat_abs / abs(y_eksak) * 100;  % galat relatif (%)

%Tampilkan hasil
fprintf('h = %.6f, n = %d\n\n', h, n);
fprintf('Tabel nilai y pada tiap langkah:\n');
fprintf('%4s %10s %12s %12s %12s %12s\n', 'i', 'x', 'Eksak', 'Euler', 'Heun', 'RK4');
for i = 1:n+1
  fprintf('%4d %10.5f %12.6f %12.6f %12.6f %12.6f\n', ...
          i-1, x(i), Ye(i), Yeu(i), Yhe(i), Yrk(i));
end

fprintf('\n=== Hasil y(2) dan galat ===\n');
fprintf('Eksak (analitik) : %.8f\n\n', y_eksak);
fprintf('%-15s %12s %14s %16s\n', 'Metode', 'y(2)', 'Galat Abstrak', 'Galat Relatif(%)');
for k = 1:3
  fprintf('%-15s %12.8f %14.8f %16.8f\n', nama{k}, hasil(k), galat_abs(k), galat_rel(k));
end

%Grafik
xf = linspace(x0, xn, 200);
figure;
plot(xf, yex(xf), 'k-', 'LineWidth', 2); hold on;
plot(x, Yeu, 'r--o', x, Yhe, 'b--s', x, Yrk, 'g--^');
grid on; xlabel('x'); ylabel('y');
legend('Eksak', 'Euler', 'Heun', 'RK4', 'Location', 'northwest');
title('Perbandingan Metode Numerik, n = 15');

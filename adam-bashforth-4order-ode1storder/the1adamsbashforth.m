%       Nama    : Z A I T U N
%       NIM     : H111 16 302
%       Tugas   : Komputasi Matematika
%       Program : Metode Adams-Bashforth Orde 4 untuk PD Orde Satu
clear;
clc;

%Persamaan Differensial Orde Satu
f = @(t,y) y - t^2 + 1;

%Solusi Analitik dari f
g = @(t) t^2 + 2*t + 1 - 0.5*exp(t);

%Diketahui
t0=0;  y0=0.5;
h=0.2;
n=10;
%==========================PERHITUNGAN==========================
for i=1:1:n+1
    t(i) = t0 + (i-1)*h;
end
%Runge Kutta 4th Order Method
y1(1) = y0;
for i=1:1:3
    k1=h*f(t(i),y1(i));
    k2=h*f(t(i)+0.5*h,y1(i)+0.5*k1);
    k3=h*f(t(i)+0.5*h,y1(i)+0.5*k2);
    k4=h*f(t(i)+h,y1(i)+k3);
    y1(i+1)=y1(i)+(1/6)*(k1+2*k2+2*k3+k4);
end

%predictor-corector
for i=4:1:n
    %Predictor
    y2p(i+1)=y1(i) + (h/24)*(55*f(t(i),y1(i))-59*f(t(i-1),y1(i-1))+37*f(t(i-2),y1(i-2))-9*f(t(i-3),y1(i-3)));
    %Corector
    y1(i+1)= y1(i) + (h/24)*(9*f(t(i+1),y2p(i+1))+19*f(t(i),y1(i))-5*f(t(i-1),y1(i-1))+f(t(i-2),y1(i-2)));
end

%Solusi Analitik
for i=1:n+1
    y4(i)=g((i-1)*h);
end

%Tabel Hasil Perhitungan
fprintf('  Table 5.14 Page 313, Adams_Bashforth 4th Pred_Cor\n');
fprintf('  ------------------------------------------\n');
     fprintf('    xi      yi         y(xi)     |y(xi)-yi|  \n');
fprintf('  ------------------------------------------\n');
for i=1:n+1
disp(sprintf('%6.1f   %8.7f   %8.7f   %8.7f%g',(i-1)*h,y1(i),y4(i),abs(y1(i)-y4(i))));
end
fprintf('  ------------------------------------------\n');

%Ploting Hasil Perhitungan
hold on
title('Graph of Adams-Bashforth 4th Method');
plot(t(1:n+1),y4(1:n+1),'g-','LineWidth',4);
plot(t(1:n+1),y1(1:n+1),'r--','LineWidth',2);

xlabel('t');ylabel('f');
legend('Analitik','Adams-Bashforth',fontsize=12);
set(legend,'location','west');
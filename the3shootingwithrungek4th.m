%              Nama    : Z A I T U N
%              NIM     : H111 16 302
%              Tugas   : Komputasi Matematika
%              Program : Metode Shooting Penyelesaian PD Orde 2
clear;
clc;

%PD Orde 2      y'' = -(2y'/x) + (2y/x^2) + sin(ln|x|)/x^2
f = @(x,y,z) z;
g = @(x,y,z) (-2*z)/x + (2*y)/x^2 + sin(log(x))/x^2;
g1 = @(x,y,z) (-2*z)/x + (2*y)/x^2;

%Solusi Analitik
c2=(1/70)*(8-12*sin(log(2))-4*cos(log(2)));
c1= (11/10)-c2;
p = @(x) c1*x + c2/x^2 - (3/10)*sin(log(x)) - (1/10)*cos(log(x));
%
x0 = 1; y0 = 1;
%
xn = 2; yn = 2;
%Tebakan
z11 = 0;
z22 = 1;
%Nilai Grid
h=0.1;
n=(xn-x0)/h;
%===============PERHITUNGAN Runge Kutta 4th Order Method===============
for i=1:1:n+1
    x(i) = x0 + (i-1)*h;
end
%Tebakan Pertama
y1(1) = y0;
z1(1) = z11;
for i=1:n
    k1y=h*f(x(i),y1(i),z1(i));
    k1z=h*g(x(i),y1(i),z1(i));
    k2y=h*f(x(i)+0.5*h,y1(i)+0.5*k1y,z1(i)+0.5*k1z);
    k2z=h*g(x(i)+0.5*h,y1(i)+0.5*k1y,z1(i)+0.5*k1z);
    k3y=h*f(x(i)+0.5*h,y1(i)+0.5*k2y,z1(i)+0.5*k2z);
    k3z=h*g(x(i)+0.5*h,y1(i)+0.5*k2y,z1(i)+0.5*k2z);
    k4y=h*f(x(i)+h,y1(i)+k3y,z1(i)+k3z);
    k4z=h*g(x(i)+h,y1(i)+k3y,z1(i)+k3z);
    y1(i+1)=y1(i)+(1/6)*(k1y+2*k2y+2*k3y+k4y);
    z1(i+1)=z1(i)+(1/6)*(k1z+2*k2z+2*k3z+k4z);
end
% Tebakan Kedua
y2(1) = 0;
z2(1) = z22;
for i=1:n
    k1y=h*f(x(i),y2(i),z2(i));
    k1z=h*g1(x(i),y2(i),z2(i));
    k2y=h*f(x(i)+0.5*h,y2(i)+0.5*k1y,z2(i)+0.5*k1z);
    k2z=h*g1(x(i)+0.5*h,y2(i)+0.5*k1y,z2(i)+0.5*k1z);
    k3y=h*f(x(i)+0.5*h,y2(i)+0.5*k2y,z2(i)+0.5*k2z);
    k3z=h*g1(x(i)+0.5*h,y2(i)+0.5*k2y,z2(i)+0.5*k2z);
    k4y=h*f(x(i)+h,y2(i)+k3y,z2(i)+k3z);
    k4z=h*g1(x(i)+h,y2(i)+k3y,z2(i)+k3z);
    y2(i+1)=y2(i)+(1/6)*(k1y+2*k2y+2*k3y+k4y);
    z2(i+1)=z2(i)+(1/6)*(k1z+2*k2z+2*k3z+k4z);
end
%Pendekatan Linier Shooting 
y3(1)=y0;
for i=2:n+1
    y3(i)= y1(i)+ (y2(i)*(yn-y1(n+1)))/y2(n+1);
end
%Solusi Analitik
for i=1:n+1
    y4(i)=p((i-1)*h+x0);
end
%Tabel Hasil Perhitungan
fprintf('  Table 11.1 Page 676, Shooting Method With Runge Kutta 4th Order\n');
fprintf('  ----------------------------------------------------------------------\n');
     fprintf('   xi       y1i          y2i       y3i~y(xi)      y(xi)     |y(xi)-y3i|\n');
fprintf('  ----------------------------------------------------------------------\n');
for i=1:n+1
disp(sprintf('%6.1f   %8.8f   %8.8f   %8.8f   %8.8f  %8.9f%g',x0+(i-1)*h,y1(i),y2(i),y3(i),y4(i),abs(y4(i)-y3(i))));
end
fprintf('  ----------------------------------------------------------------------\n');
% Plot tiga pendekatan
hold on
title('Numerical Solution with Shooting Method');
plot(x(1:n+1),y4(1:n+1),'y-','LineWidth',4);
plot(x(1:n+1),y1(1:n+1),'b--','LineWidth',2);
plot(x(1:n+1),y2(1:n+1),'g--','LineWidth',2);
plot(x(1:n+1),y3(1:n+1),'m--','LineWidth',2);
xlabel('x');ylabel('y');
legend('Solusi Analitik','Tebakan 1','Tebakan 2','Solusi pendekatan terbaik',fontsize=12);
set(legend,'location','NorthWest');
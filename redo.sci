clc;
clear;
clf;
ht = (6.626 * 10^-34) / (2 * %pi); //kg.m^2.s^-1 or Js
ht2 = ht*ht
h = 6.626 * 10^-34//
e = 1.6 * (10^(-19))//c
c = 3*10^8//m
k = 516//Nm^-1
xe = .0174//dimensionless
mu = 1.627*10^-27//kg
rin = .1274  //nm
req = rin*(10^-9)//m
//rin = input("minimum r")
rmax = req *2
rmin = rmax*.25
n = 1000
r= linspace(rmin ,rmax, n);
vm = 1/(2*xe)-.5//dimensionless
//wew = 
we  = sqrt(k/mu)//kg^-1
eps = we/(4*xe)//kg^-1
disp("eps ",eps)
dish = ht*we/(4*xe)//J????????
die = dish/e
av = sqrt(k/(2*dish))//m^-1
for i=1:4
    ev(i) = ((i-.5)*we*h*c - ((i-.5)^2)*we*xe*h*c)
    disp("Ev= ",ev(i)/e)
end
for i = 1:n
    diffa(i) = req-r(i)
    rdiffa(i)= av*diffa(i)
    inverse(i) = exp(rdiffa(i))
    v(i) = dish*((1-inverse(i))^2)
end
plot(r,v)
xlabel('x (m)');
ylabel('V(x) (J)');
title('Morse Potential')
dr = (r(n) - r(1)) / (n - 1);
a = diag(ones(1, n-2) * -2) + diag(ones(1, n-3), 1) + diag(ones(1, n-3), -1);
c = (-ht^2) / (2 * mu * dr^2); 
// Kinetic energy
T = c * a;
V = diag(v(2:n-1));
// Hamiltonian
H = T + V;
[eigenf,b] = spec(H)
count  = 0
for i=1:n-2
    if b(i,i)<dish
        count = count +1
    end
end
disp("Number of bound states: ",count)

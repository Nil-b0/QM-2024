// Task List 3: 1 D Harmonic Oscillator Program 1 and 2
clc;
clear;
clf;
m = 1.67 * 10^-27; 
w = 5.34 * 10^21;  
n = 1001;    
x = linspace((-15*10^-15), (15*10^-15), n); 
ht = (6.626 * 10^-34) / (2 * %pi); 
ht2 = ht*ht
m2 = 2*m
e = 1.6 * (10^(-19))
//potential energy
for i = 1:n
    v(i) = 0.5 * m * w^2 * x(i)^2;
end

// Plot potential energy
plot(x, v);
xlabel('x (m)');
ylabel('V(x) (J)');
title('Harmonic Oscillator Potential');
for i = 1:4
    en(i) = (i - 0.5) * ht * w; // Energy levels
end
disp("The Theratical Energy Values are: ",en)
dx = (x(n) - x(1)) / (n - 1);
a = diag(ones(1, n-2) * -2) + diag(ones(1, n-3), 1) + diag(ones(1, n-3), -1);
c = (-ht^2) / (2 * m * dx^2); 
// Kinetic energy
T = c * a;
V = diag(v(2:n-1));
// Hamiltonian
H = T + V;
[a,b] = spec(H)

xne = x(1,2:n-1)
xn = x'
xnew = xn(2:n-1,1)

vn = v(2:n-1,1)
// Calculating turning points
for i = 1 
    xt(i) = sqrt((2 * i - 1) * ht / (m * w));
end
disp("Turning Points:");
disp(xt);
a0  = (m*w/(%pi*ht))^.25
count = 4
//getting eigen values
for i = 1:count
    psi(:,i) = a(:,i)
    eigv(i) = b(i,i)
    // squaure of psi
    psi2 = psi .* psi
    //plotting the bound states and straight lines
    //plot(x, b(i)*ones(1, n), 'r')
end
disp("The Energy Eigen Values are:",eigv)
//to calculate xmid and xmid2
for i = 1:n-3
    xmid(i,1)= (xnew(i)+xnew(i+1))*.5
end
xmid2 = xnew(2:n-3)
for i = 1:count
   Nconst(i) =  inttrap(xnew(:,1),psi2(:,i))
   wf(:,i) = psi(:,i)/sqrt(Nconst(i))
   //expectation value of x
   wf2(:,i) = wf(:,i).*wf(:,i)
   yn1(:,i) = wf2(:,i).*xnew(:,1)
   xexpect(i) = inttrap(xnew(:,1),yn1(:,i))
   //expectation value of x^2
   yn2(:,i) = yn1(:,i).*xnew(:,1)
   x2expect(i) = inttrap(xnew(:,1),yn2(:,i))
  //potential expectation
   yn3(:,i) = wf2(:,i).*vn(:,1)
   vexpect(i) = inttrap(xnew(:,1),yn3(:,i))
   //momentum expectation
    scf(1)
    //window 1 is for x vs wave function
   subplot(count,1,i)
   plot(xnew,wf(:,i))
   scf(2)
   //window 2 is for x vs sqauared wave function s
   subplot(count,1,i)
   plot(xnew, wf2(:,i))
   
 end
 //momentum expectation
for j=1:count
   for i= 1:n-3
       midpsi(i,j) = (wf((i+1),j) + wf(i,j))*.5
       difpsi(i,j)  = (wf((i+1),j) - wf(i,j))/dx 
 end
 yn4(:,j) = midpsi(:,j).*difpsi(:,j)
 pexpect(j) = -(inttrap(xmid(:,1),yn4(:,j))).*%i*ht
 //p^2
for k = 1:n-4
    mid2psi(k,j) = (midpsi((k+1),j) + midpsi(k,j))*.5
    dif2psi(k,j) = (difpsi((k+1),j) - difpsi(k,j))/dx
 end
 yn5(:,j) = mid2psi(:,j).*dif2psi(:,j)
 p2expect(j) = -ht2 *(inttrap(xmid2(:,1),yn5(:,j)))
end
//to get the expectation values for kinetic and total energy to check it against given eigen values
kexpect = p2expect./(m2)
eexpect = kexpect + vexpect
disp("The Energy Expectation Values are:",eexpect)
//checking the uncertainty
varx2 = sqrt(x2expect - ((xexpect).^2))
varp2 = sqrt(p2expect -((pexpect).^2))
uncertainty = (varx2.*varp2)./(ht/2)
//to verify about large value of n
scf(3)
lpsi  =a(:,20)
lpsi2 = lpsi.*lpsi
plot(xnew,lpsi2)
z = gca()
z.data_bounds=[-15D-15,0;15D-15,0.0025]




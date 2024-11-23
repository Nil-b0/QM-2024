 //for program 1 and 2 of tasklist 4
clc
clear
clf
n= 1000//input("Enter n: ")
 m2 = (9.1*(10^(-31)))*2
e = 1.6 * (10^(-19))
e0 = 8.854*10^-12
r = linspace((10^-15 ),(10^-9),n)
dr  = ((10^-9) - (10^-15))/(n-1)
dr2 = dr^2
ht = (6.626*(10^-34))/(2*%pi)
ht2=((6.626*(10^-34))/(2*%pi))^2
c = (-ht2)/(m2*dr2*e)
for i =1:n
 v(i) = -e/(4*%pi*e0*r(i))
end
scf(0)
plot(r,v,'r')
z = gca();
z.data_bounds = [1D-15, -30; 1.2D-9, 0];
vn = v(2:n-1,1)
vd=diag(vn)
//kintectic matrix part
a = diag((ones(1,n-2)*-2))+ diag((ones(1,n-3)),1) + diag((ones(1,n-3)),-1)//identity matrix of 1xn with -2 of n*n size
//disp("The Second Derivative matrix: ")
//disp(a)
km = c*a
//disp("The Kinetic Matrix : ")
//disp(km)
h = km + vd
//disp("The Hamiltonian is: ")
//disp(h)
[A,B] = spec(h)
disp("The Grounnd state:", B(1,1))
en = -13.6
disp("The Theoratical Value for Ground State of Hydrogen is: ",en)
percent = -((abs(en-B(1,1)))*100/B(1,1))
disp("percentage error: ", percent)
show_window()
count = 4
//getting eigen values
for i = 1:count
    psi(:,i) = A(:,i)
    // squaure of psi
    psi2 = psi .* psi
end
xn = r'
xnew = xn(2:n-1,1)
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
[m,k]=max(wf2(:,1))
disp("The most probable radius for first state is: ",r(k))

 //momentum expectation
disp("The Potential Energy Expectation Values are :",vexpect)
for j=1:count
   for i= 1:n-3
       midpsi(i,j) = (wf((i+1),j) + wf(i,j))*.5
       difpsi(i,j)  = (wf((i+1),j) - wf(i,j))/dr
 end
 yn4(:,j) = midpsi(:,j).*difpsi(:,j)
 pexpect(j) = -(inttrap(xmid(:,1),yn4(:,j))).*%i*ht
 //p^2
for k = 1:n-4
    mid2psi(k,j) = (midpsi((k+1),j) + midpsi(k,j))*.5
    dif2psi(k,j) = (difpsi((k+1),j) - difpsi(k,j))/dr
 end
 yn5(:,j) = mid2psi(:,j).*dif2psi(:,j)
 p2expect(j) = -ht2 *(inttrap(xmid2(:,1),yn5(:,j)))
end
for i=1:4
  scf(3)
  subplot(4,1,i)
  rf(:,i) = psi(:,i)./xnew(:,1)
  plot(xnew,rf(:,i))
end
//to get the expectation values for kinetic and total energy to check it against given eigen values
kexpect = p2expect./(m2*e)
disp("The Kinetic Energy Expectation Values are:", kexpect)
eexpect = kexpect + vexpect
disp("The Total Energy Expectation Values are: ",eexpect)
//checking the uncertainty
varx2 = sqrt(x2expect - ((xexpect).^2))
varp2 = sqrt(p2expect -((pexpect).^2))
uncertainty = (varx2.*varp2)./(ht/2)

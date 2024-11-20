//for program 1 and 2 of tasklist 1 part A
clc
clear
n= 1000//input("Enter n: ")
x = linspace((-2*10^-10),(2*10^-10),n)
//disp(x)
v = zeros(n,1)
for i = 1:n
    if abs(x(i))<(1*(10^-10)) then
        v(i) = -100
        end
end
//disp(v)
plot(x,v,'b')
vn = v(2:n-1,1)
vd=diag(vn)
//disp("The potential energy matrix: ")
//disp(vd)
//kintectic matrix part
a = diag((ones(1,n-2)*-2))+ diag((ones(1,n-3)),1) + diag((ones(1,n-3)),-1)//identity matrix of 1xn with -2 of n*n size
//disp("The Second Derivative matrix: ")
//disp(a)
ht2=((6.626*(10^-34))/(2*%pi))^2
m2 = (9.1*(10^(-31)))*2
dx2 = ((4*(10^(-10)))/(n-1))^2//2-(-2) = 4
e = 1.6 * (10^(-19))
c = (-ht2)/(m2*dx2*e)
km = c*a
//disp("The Kinetic Matrix : ")
//disp(km)
h = km + vd
//disp("The Hamiltonian is: ")
//disp(h)
[A,B] = spec(h)
z = gca();
z.data_bounds = [-2.5*10^-10, -150; 2.5*10^-10, 50];
//get bound states
count = 0
for i = 1:n-2
    if B(i,i)<0
        count = count +1
        //to get the bound state eigen values
        eigv(i) = B(i,i) 
    end
end
disp("The nunber of bound states is : ")
disp(count)
//getting eigen values
for i = 1:count
    psi(:,i) = A(:,i)
    // squaure of psi
    psi2 = psi .* psi
    //plotting the bound states and straight lines
    plot(x, eigv(i)*ones(1, n), 'r')
end
xn = x'
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
 //momentum expectation
 dx  = (4*(10^(-10)))/(n-1)
 ht = (6.626*(10^-34))/(2*%pi)
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
kexpect = p2expect./(m2*e)
eexpect = kexpect + vexpect
//checking the uncertainty
varx2 = sqrt(x2expect - ((xexpect).^2))
varp2 = sqrt(p2expect -((pexpect).^2))
uncertainty = (varx2.*varp2)./(ht/2)


show_window()
   

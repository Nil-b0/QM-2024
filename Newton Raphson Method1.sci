clc
clear
v0 = -100*(1.6*10^-19)
a = 1*(10^-10)
m2 = (9.1*(10^(-31)))*2
ht2=((6.626*(10^-34))/(2*%pi))^2
//Value of z0
z0 = a*sqrt((m2/ht2)*abs(v0))
disp("z0 is: ")
disp(z0)
//z0 interval and root numbers
n = int(z0/(%pi/2))+1
disp("The total number of roots:")
disp(n)
//
//defiing f(z),f'(z) even solutions
function [f,df]=fdf(z)
    rz = z0/z
    sqrz= sqrt((rz)^2 -1)
    f = (tan(z) - sqrz)
    df = (sec(z)^2) + (((z0^2)/(z^3))*(1/sqrz))
endfunction
function [g,dg]=gdg(z)
    rz = z0/z
    sqrz= sqrt((rz)^2 -1)
    g = -cotg(z) - sqrz
    dg =  ((csc(z))^2) + (((z0^2)/(z^3))*(1/sqrz))
    endfunction
//defining g(z) and g'(z) odd solutions
//energy gunction
function [en] = ene(z)
en = (abs(v0)-((((z.^2)*ht2)/(m2*a*a))))/(1.6*(10^-19))   
endfunction
//ep = input("Enter accuracy for Newton Raphson:")
//zin = linspace(0,z0,(%pi/2))
z = []
en=[]
sz = size(z)
disp("The size of z")
disp(sz)
j=1
while sz <n
    //To decide on odd or even solution 
zin = input("pleae input a guess") 
      zi=[]
      zi(1)= zin
      en1(1) =  ene(zin)
      disp("look at en1:")
      disp(en1(1))
      en2 = 0 //initialization
      ck = int(zin/(%pi/2))+1
disp("ck")
disp(ck)
modu = modulo(ck,2)
disp("modulo")
disp(modu)
     for i=1:7
     if modu==1 then
     // (en1-en2)>.001
      [f,df] = fdf(zi(i))
     disp("fz and dfz:")
     disp(f)
     disp(df)
     ratio = f/df
     disp("ratio")
     disp(ratio)
     zi(i+1) = (zi(i) - (f/df))
     disp(zi(i+1))
     en1(i+1) = ene(zi(i+1))
     disp("ene2")
     en2 = en1(i+1)
     disp(en2)
 else
      //disp("Odd Solution ")
      [g,dg] = gdg(zi(i))
     //disp("gz and dgz:")
     //disp(g)
     //disp(dg)
     ratio = g/dg
     //disp("ratio")
     //disp(ratio)
     zi(i+1) = (zi(i) - (ratio))
     //disp(zi(i+1))
     en1(i+1) = ene(zi(i+1))
     //disp("ene2")
     en2 = en1(i+1)
     //disp(en2)
end
end
disp("i")
disp(i)
 disp("the root and energy")
 disp(zi(i))
 disp(en1(i))
 z(j) = zi(i)
 j= j+1  
 sz = sz+1
 en(j) = en1(i)
end 
disp("All roots:")
disp(z)
disp("all energy:")
disp(en)


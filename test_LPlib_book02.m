clear all;
close all;
clc;
A=[-2 -1 -4 -2; 5 3 1 2; 5 3 0 1;4 6 2 5];
c=[3; -4; 5; -2]; 
b=[4; 18; -13; 10];
Eqin=[-1;-1;1;0];
MinMaxLP=1;
c0=0;
% all variables is already >=0. 
% section 2.3.5 is about Transformations of Constraint Ranges and Variable Constraints  
[A, c, b, Eqin, MinMaxLP, c0] = general2canonical(A, c, b, Eqin, MinMaxLP, c0)

clear all;
close all;
clc;
A=[-2 -1 -4 -2; 5 3 1 2; 5 3 0 1;4 6 2 5];
c=[3; -4; 5; -2]; 
b=[4; 18; -13; -10];
Eqin=[-1;-1;1;1];
MinMaxLP=1;
c0=0;
% all variables is already >=0.
% section 2.3.5 is about Transformations of Constraint Ranges and Variable Constraints  
[A, c, b, Eqin, MinMaxLP] = general2standard(A, c, b, Eqin, MinMaxLP)


clear all;
close all;
clc;

A=[-2 3 5; -1 -2 2];
b=[10;-5];
c=[3; -2; 2];
Eqin=[1;1];
MinMaxLP=-1;
VarConst=[2;2;2];
[DA, Dc, Db, DEqin, DMinMaxLP, DVarConst] = primal2dual(A, c, b, Eqin, MinMaxLP, VarConst);
[A, c, b, Eqin, MinMaxLP, VarConst] = primal2dual(DA, Dc, Db, DEqin, DMinMaxLP, DVarConst);


clear all;
close all;
clc;

A=[-3 4 -3; 1 -2 4; 2 -1 -2];
b=[-3;4;6];
c=[2; 4; 0];
Eqin=[0;1;-1];
MinMaxLP=1;
VarConst=[0;2;1];
[DA, Dc, Db, DEqin, DMinMaxLP, DVarConst] = primal2dual(A, c, b, Eqin, MinMaxLP, VarConst);
% can't do that
% [A, c, b, Eqin, MinMaxLP, VarConst] = primal2dual(DA, Dc, Db, DEqin, DMinMaxLP, DVarConst)

clear all;
close all;
clc;

A=[2 0 2 3; 0 -2 -2 -6];
b=[10;-6];
c=[1; 0; -1; -3];
Eqin=[0;0];
[xsol, fval, exitflag, iterations] = rsa(A, c, b, Eqin)




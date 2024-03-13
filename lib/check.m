function boolean = check(u,p_limits)
A=p_limits*pi/180;
if(abs(u(1))<A && abs(u(2))<A && abs(u(3))<A && abs(u(4))<A)
    boolean=1;
else
    boolean=0;
end
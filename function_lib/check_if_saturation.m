function boolean = check_if_saturation(u,plim)



if(any(u<=plim(:,1)) || any(u>=plim(:,2)))
    boolean=1;
else
    boolean=0;
end
% input is quaternion q=[w,x,y,z], suppport vector
function roll=quat_to_roll(q)
w=q(:,1);
x=q(:,2);
y=q(:,3);
z=q(:,4);
dcm21 = 2 * (w .* x + y .* z);
dcm22 = w.*w - x.*x - y.*y + z.*z;
roll = atan2(dcm21, dcm22);
end
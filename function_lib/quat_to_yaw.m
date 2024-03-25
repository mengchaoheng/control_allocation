% input is quaternion q=[w,x,y,z], suppport vector
function yaw=quat_to_yaw(q)
w=q(:,1);
x=q(:,2);
y=q(:,3);
z=q(:,4);
dcm10 = 2 * (x .* y + w .* z);
dcm00 = w.*w + x.*x - y.*y - z.*z;
yaw = atan2(dcm10, dcm00);
end
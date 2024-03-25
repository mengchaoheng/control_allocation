% input is quaternion q=[w,x,y,z], suppport vector
function pitch=quat_to_pitch(q)
w=q(:,1);
x=q(:,2);
y=q(:,3);
z=q(:,4);
dcm20 = 2 * (x .* z - w .* y);
pitch = asin(-dcm20);
end
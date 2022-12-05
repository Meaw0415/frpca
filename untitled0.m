a =[1,2,3;4,5,6;8,8,8;9,9,9];
b = [3,4,5;6,7,7;4,7,3;2,8,2];

[U1,S1,z] = svd(a);
[U2,S2,z0] = svd(b);
% [j,k,l] = frpca_subspace_merge2(x,y,x0,y0);






% a=U1';
% Z = a*U2;
% [Q,R] = qr(U2-U1*Z);
% 
% temp = [S1,Z*S2;0,R*S2] ;
% [U0,Sm,Vm] = svd(temp);
% Um = [U1,Q]*U0;

//SOLUTION OF SYSTEMS OF LINEAR EQUATIONS
clc;
A = [0.531 0.196 0.236 5.032; 
0.201 0.372 2.987 0.421;  
0.217 4.691 0.279 0.237; 
3.910 0.129 0.283 0.107];
b = [0.395 0.432 0.127 0.458];
//by corollary 2 from diagonal dominance we have convergence
P = [4 3 2 1];
A = A(P,:);
b = b(P);

disp('Matrix A:');
disp(A);
disp('Vector b:');
disp(b);

eps0=0.001; //absolute error of system members
eps1=0.001; //required solution accuracy
condA=cond(A);
printf("Condition number %f\n",condA);

//compute matrix and vector norms
normA=norm(A,'fro'); //Frobenius norm
normB=norm(b,'fro');
printf("Norm A = %f\n", normA);
printf("Norm b = %f\n", normB);

absDelta=norm(inv(A))*eps0;
relDelta = condA*eps0/norm(b);
fullRelDelta = (condA/(1-(condA*(eps0/normA))))*(eps0/normA + eps0/normB);
printf("Absolute error: %f\n",absDelta);
printf("Relative error: %f,\nFull relative error: %f\n",relDelta,fullRelDelta);

if(relDelta>eps1)
 printf("Couldn't be solved within eps=%f\n",eps1);
end;

//Gauss method
printf("\nThe Gauss method");
Ab = [A b'];
disp('Extended matrix Ab:');
disp(Ab);
C=rref(Ab);
disp('Triangle matrix C:');
disp(C);
[n,m]=size(C); 
x=C(:,m);
disp('Vector X:');
disp(x);

//Simple iteration method with accuracy 0.01
printf("\nThe simple iteration method");
[M,N]=size(A);
eps1=0.01;
//transform the system to a form convenient for the method
E=eye(A); //identity matrix of size A
Diag=diag(A); //vector with diagonal elements of A
for i=1:M
  CC(i, :) = A(i, :) / Diag(i);
  B(i) = b(i) / Diag(i);
end;
C=(E-CC); //new matrix for calculations
//check convergence conditions of the iterative process
normC=norm(C,'fro');
printf("\nMatrix C norm is: %f\n",normC);
if(normC>1)
  disp(normC,'> 1, couldn't be solved using iteration')
end;
normBB=norm(B,'fro');
k=ceil(log(eps1 * (1- normC) / normBB) / log(normC)); //required number of steps for the method
printf("Steps to do: %d\n",k);
X=zeros(N,1);
for i=1:k
 XK=X;
 disp(i);
 X=C*X+B 
 delta=norm(XK-X,'fro'); //difference between steps
end;
deltaX=normC/(1 - normC) * delta; //a priori error
printf("Solution error is %f\n ",deltaX);

//Jacobi method with accuracy 0.001
printf("\nThe Jacobi method");
eps1=0.001;
for i=1:size(A,1)
  D(i,:) = E(i,:)*Diag(i);
end;
X=b'./Diag;
k=100;
for i=1:k
 XK=X;
 disp(i);
 X=(b'-(A - D)*X)./Diag 
 delta=norm(XK-X,'fro');
 printf("Delta = %f\n", delta);
 if delta < eps1 then
  break;
 end
end;
C = (A - D)/D;
normC = norm(C, 'fro');
apostDelta = (normC * delta)./(1-normC);
printf("A posteriori delta = %f\n", apostDelta);

//Seidel method with accuracy 0.001
printf("\nThe Seidel method");
L=tril(A, -1);
R=triu(A, 1);
invD=inv(D);
k=100;
for i=1:k
 XK=X;
 disp(i);
 X = inv(E + invD*L)*invD*(b' - R*X)
 delta=norm(XK-X,'fro');
 printf("Delta = %f\n", delta);
 if delta < eps1 then
  break;
 end
end;
C = inv(E + invD*L)*invD*R;
normC = norm(C, 'fro');
apostDelta = (normC * delta)./(1-normC);
printf("A posteriori delta = %f\n", apostDelta);

syms y x 


rho=1; sigma=1; omega=1; mu=4*pi*10^(-7); 

% xnew=2; ynew=0; 
% z=1;
syms xnew ynew 
syms z

% definitions
R=rho*sin(y);
zRelative=z-rho*cos(y);
k=(mu*sigma*omega*R^2/(4*pi));

posthetta=asin(sqrt(xnew.^2+ynew.^2)/sqrt(xnew.^2+ynew.^2+z^.2));
posphi=asin(ynew/sqrt(xnew.^2+ynew.^2));

% biot savart variables
BsphereInti=  k*(zRelative)* ...
    (R*cos(x))/(xnew^2+ynew^2+(zRelative)^2+ ...
    (R)^2-2*R*(ynew*sin(x)+xnew*cos(x)))^(3/2);

BsphereIntj= k*(zRelative)* ...
    (R*sin(x))/(xnew^2+ynew^2+(zRelative)^2+ ...
    (R)^2-2*R*(ynew*sin(x)+xnew*cos(x)))^(3/2);

BsphereIntk= k*(R^2-R*(ynew*sin(x)+xnew*cos(x)))/(xnew^2+ynew^2+(zRelative)^2+ ...
    (R)^2-2*R*(ynew*sin(x)+xnew*cos(x)))^(3/2);

% Vector potential variables
Br = (((mu * rho.^4 * omega * sigma)/3) * (1/(sqrt(z^2+xnew.^2+ynew.^2)).^3)) * 2*cos(posthetta);
Btheta = (((mu * rho^4 * omega * sigma)/3) * (1/(sqrt(z^2+xnew.^2+ynew.^2)).^3)) * sin(posthetta);

BxVector=Br.*sin(posthetta).*cos(posphi) + Btheta.*cos(posthetta).*cos(posphi);
ByVector= Br.*sin(posthetta).*sin(posphi) + Btheta.*cos(posthetta).*sin(posphi);
BzVector = Br.*cos(posthetta) - Btheta.*sin(posthetta);

%Storage arrays
VectorArray=zeros(3,100);
BiotSavartArray=zeros(3,100);

for i = 1:10
%integration of biot savart

iComponent = matlabFunction(subs(subs(subs(BsphereInti,ynew,1),xnew,i),z,1));
jComponent = matlabFunction(subs(subs(subs(BsphereIntj,ynew,1),xnew,i),z,1));
kComponent = matlabFunction(subs(subs(subs(BsphereIntk,ynew,1),xnew,i),z,1));

BiotSavartArray(1,i)=integral2(iComponent,0,2*pi,0,pi);
BiotSavartArray(2,i)=integral2(jComponent,0,2*pi,0,pi);
BiotSavartArray(3,i)=integral2(kComponent,0,2*pi,0,pi);


% Analytical solution using vector potential

VectorArray(1,i) = subs(subs(subs(BxVector,ynew,1),xnew,i),z,1);
VectorArray(2,i) = subs(subs(subs(ByVector,ynew,1),xnew,i),z,1);
VectorArray(3,i) = subs(subs(subs(BzVector,ynew,1),xnew,i),z,1);
end

iteration=1:1:10;

figure(1)
plot(iteration,VectorArray(1,iteration))
hold on
plot(iteration,BiotSavartArray(1,iteration))

legend('Vector Potential','Biot Savart')
title('X component')
xlabel('steps') 
ylabel('Magnetic Field Intensity') 

figure(2)
plot(iteration,VectorArray(2,iteration))
hold on
plot(iteration,BiotSavartArray(2,iteration))
legend('Vector Potential','Biot Savart')
title('Y component')
xlabel('steps') 
ylabel('Magnetic Field Intensity') 

figure(3)
plot(iteration,VectorArray(3,iteration))
hold on
plot(iteration,BiotSavartArray(3,iteration))
legend('Vector Potential','Biot Savart')
title('Z component')
xlabel('steps') 
ylabel('Magnetic Field Intensity') 


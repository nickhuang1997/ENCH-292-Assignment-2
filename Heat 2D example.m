% Example MATLAB script to solve the 2D heat equation using a finite
% difference discretisation

% Equation to solve is:

% d^2T/dx^2 + d^2T/dy^2 = 0

% The numerical solution is obtained by solving the discretised form of
% this equation:

% (T_m+1,n - 2*T_m,n + T_m-1,n)/Dx^2 + (T_m,n+1 - 2*T_m,n + T_m,n-1)/Dy = 0

% We will develop the solution for constant boundary conditions and a 3 by
% 3 set of temperatures.
% Let the temperatures be given by:
% T(1,1) T(2,1) T(3,1)
% T(1,2) T(2,2) T(3,2)
% T(1,3) T(2,3) T(3,3)
% by symmetry T(1,1) = T(3,1); T(3,2) = T(1,2); T(3,3) = T(1,3)
% Let T(1,1) = T1; T(1,2) = T2; T(1,3) = T3;
% T(2,1) = T4; T(2,2) = T5; T(2,3) = T6
% We can then define our matrix equations:
A = zeros(6,6); % 6 unknown temperatures, and 6 equations to define these.

A(1,:) = [4 -1 0 -1 0 0];
A(2,:) = [-1 4 -1 0 -1 0];
A(3,:) = [0 -1 4 0 0 -1];
A(4,:) = [-2 0 0 4 -1 0];
A(5,:) = [0 -2 0 -1 4 -1];
A(6,:) = [0 0 -2 0 -1 4];

b=zeros(6,1);
b(1) = 0 + 100;
b(2) = 0;
b(3) = 0 + 0;
b(4) = 100;
b(5) = 0;
b(6) = 0;

% Solve the equations
T = A\b;                

% reshape temperature to generate 2D map
T = reshape(T,[3 2]);
% replicate first column from symmetry
T(:,3)=T(:,1);

% add BC cells
TpBC = zeros(5,5);
TpBC(1,:)=100;
TpBC(2:4,2:4)=T;
TpBC(1,1)=50;
TpBC(1,5)=50;

imagesc(TpBC,[0 110])
c = colorbar;
c.Label.String = 'Temperature (?C)';
colormap hot
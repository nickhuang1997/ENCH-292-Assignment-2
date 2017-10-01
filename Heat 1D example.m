% Example MATLAB script to solve the 1D heat equation using a finite
% difference discretisation

% The equation to solve is:

% d^2T/dx^2 = 0

% The numerical solution is obtained by solving the discretised form of
% this equation:

% T_i+1 - 2*T_i + T_i-1 = 0

% We will develop the solution for three cases:
% (1) constant boundary conditions
% (2) constant heat flux BC at one end and constant temperature BC at the other end.
% (3) convective heat flux BC at one end and constant temperature BC at the other end.


N = 50; % number of final cell (so N + 1 temperatures in total)
TL = 100; % ?C Constant temperature BC at hot end
TR = 20; % ?C constant temperature BC at cold end
qN = 5000; %W/m2 % heat flux at cold end
h = 10; % W /m2 /K heat transfer coefficient for convective BC
Te = 20;    % ?C ambient temperature for convective BC

k = 12; %W/m/K % thermal conductivity
L = 10; %cm % length of spoon in cm

%% simple implementation with constant temperature BCs
% this implementation solves for all N + 1 temperatures.

% Set up matrix to describe system of discretised equations
A = zeros(N+1,N+1);       % initialise A matrix with zeros
x = linspace(0,L,N+1);    % set up horizontal coordinates for each temperature

A(1,1) = 1;               % Set first row of A matrix for constant temperature BC
A(N+1,N+1) = 1;           % Set last row of A matrix for constant temperature BC

for iN = 2:N
    A(iN,iN-1:iN+1)=[1 -2 1]; % set values of A in each of the central rows
end

% Set up vector for RHS of system of discretised equations
b = zeros(N+1,1);       % set b
b(1)=TL;                % set constant temperature BC for hot end
b(N+1)=TR;              % set constant temperature BC for cold end

% Solve the equations
T = A\b;                

% Plot result
plot(x,T,'-','LineWidth', 2)        % set linewidth to ensure visible      
xlabel('Position (cm)');            % set horizontal axis label
ylabel('Temperature (?C)');         % set vertical axis label
axis([0 10 0 100]);                 % define axis limits
hold all;                           % ensure subsequent plots appear in same figure
%% simple implementation with a constant flux BC

% Set up matrix to describe system of discretised equations
A = zeros(N+1,N+1);     % initialise A matrix with zeros
x=linspace(0,L,N+1);    % set up horizontal coordinates

A(1,1) = 1;             % Set first row of A matrix
A(N+1,N:N+1) = [-1; 1];  % set final row of A matrix

for iN = 2:N
    A(iN,iN-1:iN+1)=[1 -2 1]; % set values of A in each of the central rows
end

% Set up vector for RHS of system of discretised equations
b = zeros(N+1,1);       
b(1)=TL;                        % set constant temperature BC for hot end
b(N+1)=-q/k*(x(2)-x(1))/100;    % set constant heat flux BC for cold end
                                % Note: factor of 100 is to convert Dx from cm to m

% Solve the equations
T = A\b;                

% Plot result
plot(x,T,'-','LineWidth', 2)               
xlabel('Position (cm)');     
ylabel('Temperature (?C)'); 
axis([0 10 0 100]);
%% simple implementation with a single convective flux BC

% Set up matrix to describe system of discretised equations
A = zeros(N+1,N+1);     % initialise A matrix with zeros
x=linspace(0,L,N+1);    % set up horizontal coordinates
Bi = h*(x(2)-x(1))/k;   % Biot number for grid.

A(1,1) = 1;             % Set first row of A matrix
A(N+1,N:N+1) = [1 -(1+Bi)];  % set final row of A matrix

for iN = 2:N
    A(iN,iN-1:iN+1)=[1 -2 1]; % set values of A in each of the central rows
end

% Set up vector for RHS of system of discretised equations
b = zeros(N+1,1);       
b(1)=TL;                        % set constant temperature BC for hot end
%b(N+1)=-q/k*(x(2)-x(1))/100;    % set constant heat flux BC for cold end
b(N+1)=-Bi*Te;    % set constant heat flux BC for cold end
                                % Note: factor of 100 is to convert Dx from cm to m

% Solve the equations
T = A\b;                

% Plot result
plot(x,T,'-','LineWidth', 2)               
xlabel('Position (cm)');     
ylabel('Temperature (?C)'); 
axis([0 10 0 100]);

% compare to analytical result:
U = (L/k + 1/h)^-1;
q = U*(TL-Te);
TL = q/h+Te

% alternatively, check the heat flux:
q_num = -k*(T(end) - T(end-1))/(x(end) - x(end - 1));


%% more efficient implementation
% 

% Set up matrix to describe system of discretised equations
A = zeros(N-1,N-1);     % initialise A matrix with zeros
                        % Note: only N-1 elements as we are going to build
                        % the constant temperature BCs into the b matrix
                        % to reduce the number of unknowns
x=linspace(0,L,N+1);    % set up horizontal coordinates

A(1,1:2) = [-2 +1];     % Set first row of A matrix
A(N-1,N-2:N-1) = [+1 -2];

for iN = 2:N-2;
    A(iN,iN-1:iN+1)=[1 -2 1]; % set values of A in each of the central rows
end

% Set up vector for RHS of system of discretised equations
b = zeros(N-1,1);       
b(1)=-T0;               % set constant temperature BC for hot end. 
                        % Note this is now set for the cell 1 so the BC 
                        % has changed, as has the corresponding row of A. 
b(N-1)=-TN;             % set constant temperature BC for cold end, as for hot end.

% Solve the equations
T = A\b;                

% Plot result
plot(x(2:end-1),T,'-','LineWidth', 2)   % Note we only have solutions for N-1 cells,
xlabel('Position (cm)');                % hence the changed coordinates of plot command
ylabel('Temperature (?C)');             
axis([0 10 0 100]);
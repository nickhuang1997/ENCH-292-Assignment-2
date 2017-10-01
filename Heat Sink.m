% Example MATLAB script to solve the 2D heat equation for a computer heat sink 
%%
% R=10;C=15
% 
% C =
% 
%     15
% 
% Adum = zeros(R,C); % Dummy 2D matrix containing weighting at each
% imagesc(Adum)
% Adum(6,8)=-4;Adum(5,8)=1;Adum(7,8)=1;Adum(6,9)=1;Adum(6,7)=1;
% imagesc(Adum)
% plot(Adum(:))
% Adum = zeros(R,C); % Dummy 2D matrix containing weighting at each cell

%%
% using a finite difference discretisation
%
% Assumes symmetric fin with spacing at base of L, thickness of base h, 
% thickness of the fin W and height of fin H
%
%       _____           
%   ^        |          
%   |        |    q=h(Ts-Tair) 
%   |   <--->|          
%   |    W/2 |          
% H |        ______________         
%   |                        ^
%   |  <------------------   | h
%   |             L/2        |
%   |   ___________________  |  
%
%       q_chip

% Equation to solve is:

% d^2T/dx^2 + d^2T/dy^2 = 0

% The numerical solution is obtained by solving the discretised form of
% this equation on a uniform grid with Dx = Dy:

% T_m+1,n - 4*T_m,n + T_m-1,n + T_m,n+1 + T_m,n-1 = 0

% The boundary conditions are defined flux conditions:
% for external surfaces q = hc(Ts - Tair)
% for base in contact with chip q = 15000 W m^-2
%%%%%%
%Redacted
%%%%%%

% set boundary conditions
q_chip = 15000; % W m^-2  - Approximate heat flux from core i7 chip

% set parameters of materials
L = 0.004; %m
W = 0.0016; %m
h = 0.002; %m
H = 0.032; %m

hc = 20; % W m^-2 K^-1
Dx = 0.05e-3 ; % m
k = 200; % W m^-1 K^-1 - Thermal conductivity of common aluminium alloy
Tair = 20; % degC
Bi = hc*Dx/k;

% set number of cells in each direction:
R = int16(H/Dx); %int16 forces integer values. There are better ways to do this.
C = L/2/Dx;

%%
% initialise geometry

geom = zeros(R,C);
geom(1:h/Dx,1:end)=1;
geom(1:end,1:W/2/Dx)=1;
idx=find(geom==1);  %identifies all cells within the domain (i.e. excludes cells containing air)

figure(1);imagesc(geom)
title('geometry')
colorbar
set(gca,'YDir','normal')

%% Set up matrices

% initialise A matrix. One row for each unknown temperature within the
% solid. Cells containing air excluded.
A = sparse(length(idx),length(idx));
%A = zeros(length(idx),length(idx));
Adum = zeros(R,C); % Dummy 2D matrix containing weighting at each 
% row and cell index. Used to assign weights to correct cells in A matrix
% when Adum is converted to a vector form it comprises one row of the A
% matrix.

b = Adum(idx); % initialise b as matrix of zeros corresponding to geometry of plate

% set discretised equations into A matrix and b vector
for iN = 1:length(idx)
    Adum = zeros(R,C); % Dummy 2D matrix containing weighting at each cell
    Adum(idx(iN)) = 2; % set current entry to a value of 2
    [iR, iC]=find(Adum==2); % find locations of the entry identified on above row
        if iR == 1 && iC == 1 % Bottom LH corner cell
            cas = 1;
            Adum(iR,iC:iC+1) = [2 -1];
            Adum(iR+1,iC) = -1;
            b(iN) = q_chip*Dx/k;
%%%%%
% Redacted
%%%%%
        
        end
%         figure(2);imagesc(Adum) % These commands can be helpful to debug.
%         set(gca,'YDir','normal')
        
        A(iN,:)=Adum(idx); % convert Adum to vector and place as row in A matrix
     
end


%% solve
tic
T = A\b; % solves system of equations. Produces temperature as a vector T
toc
Tall=zeros(R,C);
Tall(idx)=T; % reshape T vector to correspond to physical positions.


%% display result
Tall(Tall==0)=1000;
figure(3);imagesc(Tall,[5*floor(min(T(:))/5) 5*ceil(max(T(:))/5)]) % display result
set(gca,'YDir','normal')
c = colorbar;
c.Label.String = 'Temperature (?C)';
title('Temperature')
colormap hot

%% calculate heat flux
% calculates matrix for the heat flux. Set up as a 3 dimensional matrix.
% The first two dimensions correspond to real space positions. The third
% dimension has two entries corresponding to q_x and q_y.
%%%%%%
% Redacted
%%%%%%

figure(2);
cmap = colormap('jet');
cmap(end,:)=[1 1 1];
colormap(cmap);
imagesc(temp,[-2 2]*1e3)
title('q_x')
c = colorbar;
c.Label.String = 'q_x (W m^{-2} K^{-1})';
set(gca,'YDir','normal')

% Total heat loss through fin

%%%%%%
% Redacted
%%%%%%

figure(4); close; figure(4); hold off
starts= [1:5:30 40:20:400];
ends=1:5:40;
streamline(qvec(:,:,2),qvec(:,:,1),[1:5:40],ones(size(1:5:40)));%ones(size(starts)),starts,'r')
Fig1Ax1 = get(4,'Children');
Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
set(Fig1Ax1Line1, 'LineWidth', 1);
hold all
contour(Tall,[65:0.2:67 67:0.05:68],'r','LineWidth',1)
title('stream lines for heat flow')
%ENCH292 Heat and Mass Transfer Assignment 3
%Jackson Olds and Georgia Byrne

%Solving a 2D system using numerical methods
%d^2T/dx^2 + d^2T/dy^2 = 0

%Target area is 0.1 x 0.3 m^2 section containing half of one of the wires
%and half of the area between two wires.

clc 
clear all
close all

%Boundary conditions
n = 300;
x = 0.3;
y = 0.1;

R = y*n;
C = x*n;

dx = x/n;

%Known and Required values
Te = 20; %degrees celcius
h = 20; %W m-2 K -1
L = 0.3; %m
k = 0.12; %W m-1
Bi = h*dx/k; 
powerTotal = 4000; %W
surfaceArea = 6.4; % m^2
qs = powerTotal/surfaceArea; % Wm^-2

A = zeros((R-1)*(C-1),(R-1)*(C-1));
Adum = zeros((R-1),(C-1)); % Dummy 2D matrix containing weighting at each row and cell index. Used to assign weights to correct cells in A matrix
% when Adum is converted to a vector form it comprises one row of the A
% matrix.

% set discretised equations into a matrix
for iR = 1:R-1;
    for iC = 1:C-1;
        
        Adum(:) = 0; %Initialising every cell to 0
             
        if iR == 1 && iC == 1; %Top LH corner
            Adum(iR:iR+1, iC) = [2+Bi -1];
            Adum(iR, iC+1) = -1;
            
        elseif iR == 1 && iC == C-1; %Top RH corner
            Adum(iR:iR+1, iC) = [2+Bi -1];
            Adum(iR, iC-1) = -1;
            
        elseif iR == R-1 && iC == 1; %Bottom LH corner
            Adum(iR-1:iR, iC) = [-1, 2];
            Adum(iR, iC+1) = -1;
            
        elseif iR == R-1 && iC == C-1; %Bottom RH corner
            Adum(iR-1:iR, iC) = [-1, 2];
            Adum(iR, iC-1) = -1;
                        
        elseif iR == 1; %Top face
            Adum(iR, iC-1:iC+1) = [-1/2 2+Bi -1/2];
            Adum(iR+1, iC) = -1;
        
        elseif iR == R-1; %Bottom face
            Adum(iR, iC-1:iC+1) = [-1 4 -1];
            Adum(iR-1, iC) = -2;

            
        %Area around the wire
        
        elseif iR == R/2 && iC == 1;    %Top Left Corner Wire
            Adum(iR, iC:iC+1) = [-4 1];   
            Adum(iR-1, iC) = 3;
        
        elseif iR == R/2 && iC == C/30; %Top Right Corner Wire
            Adum(iR-1:iR+1, iC) = [1 -4 1];
            Adum(iR, iC-1) = 1;
            Adum(iR, iC+1) = 1;
            
        elseif iR == 0.7*R && iC == 1;  %Bottom Left Corner Wire
            Adum(iR, iC:iC+1) = [-4 1];    
            Adum(iR+1, iC) = 3;
            
        elseif iR == 0.7*R && iC == C/30; %Bottom Right Corner Wire
            Adum(iR-1:iR+1, iC) = [1 -4 1];
            Adum(iR, iC-1) = 1;
            Adum(iR, iC+1) = 1;            
             
        elseif iR == R/2 && iC < C/30; %Top of Wire
            Adum(iR, iC-1:iC+1) = [-1 4 -1];
            Adum(iR-1, iC) = -2;
%             
        elseif iR == 0.7*R && iC < C/30; %Bottom of wire
            Adum(iR, iC-1:iC+1) = [-1 4 -1];
            Adum(iR+1, iC) = -2;
%             
        elseif R/2 < iR < 0.7*R && iC == C/30; %RHS of wire
            Adum(iR-1:iR+1, iC) = [-1 4 -1];
            Adum(iR, iC+1) = -2;
            
        %Remaining sides
        
        elseif 2 <= iR < R/2 && iC == 1; %Upper LHS face
            Adum(iR-1:iR+1, iC) = [-1 4 -1];
            Adum(iR, iC+1) = -2;
            
        elseif 0.7*R < iR < R-1 && iC == 1; %Lower LHS face
            Adum(iR-1:iR+1, iC) = [-1 4 -1];
            Adum(iR, iC+1) = -2;
            
        elseif iC == C-1; %RHS symmetrical face
            Adum(iR-1:iR+1, iC) = [-1 4 -1];
            Adum(iR, iC-1) = -2;
            
        else
            Adum(iR, iC-1:iC+1) = [1 -4 1]; %All remaining internal cells
            Adum(iR-1, iC) = 1;
            Adum(iR+1, iC) = 1;
%             
        end
        
        A((iC-1)*(R-1)+iR,:)=Adum(:); % convert Adum to vector and place as row in A matrix
    end
end

Adum(:) = 0; 
b = Adum; %Set up of 'b' matrix using dummy A matrix for reference

b(1, :) = Bi*Te; %Top side to room
b(R/2:0.7*R, C/30) = 2*qs*dx/k; %RHS of wire
b(R/2, 1:C/30) = 2*qs*dx/k; %Top of wire
b(0.7*R, 1:C/30) = 2*qs*dx/k; %Bottom of wire

b = b(:); %Reshape as a vector

T = A\b; % solves system of equations. Produces temperature as a vector T
T=reshape(T,[R-1,C-1]); % reshape T vector to correspond to physical positions.

Tmax = max(T(:)); %Finding maximum value of the matrix

T(0.5*R:0.7*R, 1:C/30) = Tmax; %Set cells inside the wire to maximum temperature

tempDiff = T(1, 1) - T(1, end); %Difference between the max and min temps 
                                %of the visible floor surface
Tmin = min(T(:)); %Minimum temperature for comparison

mirrorT = fliplr(T); %Completing both sides of the wire
newT = [mirrorT T];
floorPlan = [newT newT newT newT newT newT newT newT newT newT]; %New array for all 10 wires

hold on
figure(1)
imagesc(flipud(newT),[15 40]) % display result
c = colorbar;
c.Label.String = 'Temperature (?C)';
colormap hot
title('Single Wire Temperature Map')
%xlim([0 180])


figure(2)
imagesc(floorPlan,[15 40]) % display result
c = colorbar;
c.Label.String = 'Temperature (?C)';
colormap hot
title('Full Temperature Map of Floor Space')
%xlim([0 1800])


fprintf('The maximum temperature difference across the top of the floor is %.2f degrees.\n', tempDiff)
fprintf('The maximum temperature of the floor is %.2f degrees.\n', Tmax)
fprintf('The minimum temperature of the floor is %.2f degrees. \n', Tmin)
    
    
    
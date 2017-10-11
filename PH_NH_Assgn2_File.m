clear all
%   |\\\\\\\\\\\\\\\\\\\\\\\\\\\|  
%   |___________________________|   _
%   |                           |   |
%   |                           |   |      
%   |                           |   |      
%   |_____                      |   |  H
%   |     |                     |   |
%   |_____|                     |   |
%   |                           |   |
%   |___________________________|   _
%               L

%%
%setting parameters of materials
w = .033; %m
h = .04; %m
Dx = 0.0005;       %.5 mm cell width and height
k_c       = 385;     %thermal conductivity of copper wire   %%%%% engr toolbox
k_i       = 2;     %thermal conductivity of insulation
k_t       = 4;     %thermal conductivity of tile
k_air     = 0.024;

%400/.0002;
q_wire    = 400/0.000025;

%Bi %??? do we need this?

R = int16(h/Dx);
C = w/Dx;       
%%
% initialising the geometry

geom = zeros(R,C);
% model = zeros(R,C);

geom(1:end,1:end) = 1;

% model(5:10,1:2.5)  = 4;
% model(40:80,1:end) = 6;   

idx=find(geom==1);      %finds all the cells with 1's 
                        %(i.e. in the domain. should be all the cells in our case)
% figure(1);imagesc(geom)
% title('geometry')
% colorbar
% set(gca,'YDir','normal')

% figure(2);imagesc(model)
% title('geometry')
% colorbar
% set(gca,'YDir','normal')
%%
% Setting up Matricies
%A = sparse(length(idx),length(idx));
A = zeros(length(idx),length(idx));
Adum = zeros(R,C); % Dummy 2D matrix containing weighting at each 
% row and cell index. Used to assign weights to correct cells in A matrix
% when Adum is converted to a vector form it comprises one row of the A
% matrix.

b = Adum(idx); % initialise b as matrix of zeros corresponding to geometry of plate
% figure(2);imagesc(b)
% title('Dummy matrix')
% colorbar
% set(gca,'YDir','normal')

for iN = 1:length(idx)
    Adum = zeros(R,C);          % Dummy 2D matrix containing weighting at each cell
    Adum(idx(iN)) = 2;          % set current entry to a value of 2
    [iR, iC]=find(Adum==2);     % find locations of the entry identified on above row

%%%%%%%%%%%%%%%%%%_corners of model_%%%%%%%%%%%%%%%%%%%%

        if iR == 1 && iC == 1   % Bottom LH corner cell
            cas = 1;
            Adum(iR,iC:iC+1) = [-2 1];
            Adum(iR+1,iC)    = 1;
            %b(iN) = q_wire*Dx^2/k_i;    %done
% 
%         elseif iR == 80 && iC == 1      %top left of model - exposed to air 
%             cas = 2;
%             Adum(iR,iC:iC+1) = [-(2*k_t+k_air) k_t];
%             Adum(iR+1,iC) = k_air;
%             Adum(iR-1,iC) = k_t;
%             %b(iN) = q_wire*Dx^2;    %done - with inspection of case 2
%         
%         elseif iR == 80 && iC == 66      %top right of model
%             cas = 3;
%             Adum(iR,iC-1:iC) = [k_t -(2*k_t+k_air)];
%             Adum(iR+1,iC) = k_air;
%             Adum(iR-1,iC) = k_t;
            %b(iN) = q_wire*Dx^2;    %done - with inspection of case 2  
            
%         elseif iR == 1 && iC == 66      %bottom right of model
%             cas = 4;
%             Adum(iR,iC-1:iC) = [1 -2];
%             Adum(iR+1,iC)    = 1;
%             %b(iN) = q_wire*Dx^2/k_i;    %done by inspection of cas #1
%         
% %%%%%%%%%%%%%%%%%%_where insulation meets tile_%%%%%%%%%%%%%%%%%%%%
% %                      right and left side                        %
%         
%         elseif iR == 40 && iC == 1      %left side insulation meets tile
%             cas = 5;
%             Adum(iR,iC:iC+1) = [-(2*k_t+k_i) k_t];
%             Adum(iR+1,iC)    = k_t;
%             Adum(iR-1,iC)    = k_i;
%             %b(iN) = -q_wire*Dx^2;  %done 
% 
%         elseif iR == 40 && iC == 66      %left side insulation meets tile
%             cas = 6;
%             Adum(iR,iC-1:iC) = [k_t -(2*k_t+k_i)];
%             Adum(iR+1,iC)    = k_t;
%             Adum(iR-1,iC)    = k_i;
%             %b(iN) = -q_wire*Dx^2;  %done
%             
%         elseif iR == 40 && iC > 1 && iC < 66 %do range for bottom of tile
%             cas = 24;
%             Adum(iR,iC-1:iC+1) = [k_t -(3*k_t+k_i) k_t];
%             Adum(iR+1,iC)    = k_t;
%             Adum(iR-1,iC)    = k_i;
%             %b(iN) = -q_wire*Dx^2;  %done
%             
% %%%%%%%%%%%%%%%%%%_corners of wire _%%%%%%%%%%%%%%%%%%%%
%  
        elseif iR == 10 && iC == 1      %bottom left of wire
            cas = 7;
            Adum(iR,iC:iC+1) = [-(2*k_c+k_i) k_c];
            Adum(iR+1,iC) = k_c;
            Adum(iR-1,iC) = k_i;
            b(iN) = -q_wire*Dx^2;   %done by inspection 
        
        elseif iR == 20 && iC == 1      %top left of wire
            cas = 8;
            Adum(iR,iC:iC+1) = [-(2*k_c+k_i) k_c]
            Adum(iR+1,iC) = k_i;
            Adum(iR-1,iC) = k_c;
            b(iN) = -q_wire*Dx^2;  %done   
        
        elseif iR == 10 && iC == 5      %bottom right corner of the wire
            cas = 9;
            Adum(iR,iC-1:iC+1) = [k_c -(2*k_c+2*k_i) k_i];
            Adum(iR+1,iC)    = k_c;
            Adum(iR-1,iC)    = k_i;
            b(iN) = -q_wire*Dx^2;  %done
            
        elseif iR == 20 && iC == 5      %top right corner of the wire
            cas = 10;
            Adum(iR,iC-1:iC+1) = [k_c -(2*k_c+2*k_i) k_i];
            Adum(iR+1,iC)    = k_i;
            Adum(iR-1,iC)    = k_c;
            b(iN) = -q_wire*Dx^2;  %done
% 
% %%%%%%%%%%%%%%%%%%_top and bottom of wire_%%%%%%%%%%%%%%%%%%%%
%         
        elseif iR == 20 && iC < 6 %do range for top of copper
            cas = 11;
            Adum(iR,iC-1:iC+1) = [k_c -(3*k_c+k_i) k_c];
            Adum(iR+1,iC)    = k_i;
            Adum(iR-1,iC)    = k_c;
            b(iN) = -q_wire*Dx^2;  %done
% 
        elseif iR == 10 && iC <6    %do range for bottom of copper
            cas = 12;
            Adum(iR,iC-1:iC+1) = [k_c -(3*k_c+k_i) k_c];
            Adum(iR+1,iC)    = k_c;
            Adum(iR-1,iC)    = k_i;
            b(iN) = -q_wire*Dx^2;  %done    
% 
% %%%%%%%%%%%%%%%%%%_top of model_%%%%%%%%%%%%%%%%%%%%
%        
%         elseif iR == 80 && iC > 0   %do range for top of t exposed to air
%             cas = 13;
%             Adum(iR,iC-1:iC+1) = [k_t -(3*k_t+k_air) k_t];
%             Adum(iR+1,iC)    = k_air;
%             Adum(iR-1,iC)    = k_t;
%             %b(iN) = -q_wire*Dx^2;  %done
% 
% %%%%%%%%%%%%%%%%%%_interior cells_%%%%%%%%%%%%%%%%%%%%
%             
        elseif iR > 10 && iR < 20 && iC < 6 && iC > 1 %range for all interior wire cells 
            cas = 14;
            Adum(iR,iC-1:iC+1) = [1 -(3) 1];
            Adum(iR+1,iC)    = 1;
            Adum(iR-1,iC)    = 1;
            b(iN) = -q_wire*Dx^2/k_c;  %done
%                     
%         elseif iR > 1 && iR < 38 && iC < 65 && iC > 1%range for all interior insulation cells 
%             %use or statement
% 
% %             if iR >8 && iR < 22 && iC < 7 
% %                 break % put at end??
% %                 %inside = False %if cells are the wire... pass this statement
% %             end
%             cas = 15;
%             Adum(iR,iC-1:iC+1) = [1 -(3) 1];
%             Adum(iR+1,iC)    = 1;
%             Adum(iR-1,iC)    = 1;
%             %b(iN) = -q_wire*Dx^2/k_i;  %done         
%         
%         elseif  iR > 39 && iC > 1 && iC < 66 %range for all interior tile cells 
%             cas = 16;
%             Adum(iR,iC-1:iC+1) = [1 -(3) 1];
%             Adum(iR+1,iC)    = 1;
%             Adum(iR-1,iC)    = 1;
%             %b(iN) = -q_wire*Dx^2/k_t;  %done
% 
% %%%%%%%%%%%%%%%%%%_bottom of insulation_%%%%%%%%%%%%%%%%%%%%
% 
%         elseif iR < 2 && iC > 0 %do range for all bottom of the insulation            
%             cas = 17;
%             Adum(iR,iC-1:iC+1) = [k_i -(3*k_i) k_i];
%             Adum(iR+1,iC)    = k_i;
%             %b(iN) = -q_wire*Dx^2;  %done
% 
% %%%%%%%%%%%%%%%%%%%%%_left walls_%%%%%%%%%%%%%%%%%%%%%%%
% 
%         elseif (iR > 0 && iR < 9 && iC < 2) | (iR > 21 && iR < 39 && iC < 2)  %do range for all left walls of insulation
%             cas = 18;
%             Adum(iR,iC:iC+1) = [-3 k_i];
%             Adum(iR+1,iC)    = k_i;
%             Adum(iR-1,iC)    = k_i;
%             %b(iN) = -q_wire*Dx^2;  
% 
%         elseif iR > 39 && iR < 80 && iC < 2%do range for all left walls of tile
%             cas = 19;
%             Adum(iR,iC:iC+1) = [-3 k_t];
%             Adum(iR+1,iC)    = k_t;
%             Adum(iR-1,iC)    = k_t;
%             %b(iN) = -q_wire*Dx^2;
%             
        elseif iR > 9 && iR < 13 && iC == 1 %do range for all left walls of copper
            cas = 20;
            Adum(iR,iC:iC+1) = [-3 k_c];
            Adum(iR+1,iC)    = k_c;
            Adum(iR-1,iC)    = k_c;
            b(iN) = -q_wire*Dx^2;
%             
% %%%%%%%%%%%%%%%%%%%%%_right walls_%%%%%%%%%%%%%%%%%%%%%%%
% 
%         elseif iR > 1 && iR < 40 && iC > 65 %do range for all right walls of insulation
%             cas = 21;
%             Adum(iR,iC-1:iC) = [k_i -3];
%             Adum(iR+1,iC)    = k_i;
%             Adum(iR-1,iC)    = k_i;
%             %b(iN) = -q_wire*Dx^2;  
% 
%         elseif iR > 39 && iR < 80 && iC > 65%do range for all right walls of tile
%             cas = 22;
%             Adum(iR,iC-1:iC) = [k_t -3];
%             Adum(iR+1,iC)    = k_t;
%             Adum(iR-1,iC)    = k_t;
%             %b(iN) = -q_wire*Dx^2;
%             
        elseif iR > 9 && iR < 11 && iC < 7 && iC > 5 %do range for all right walls of copper
            cas = 23;
            Adum(iR,iC-1:iC) = [k_c -3];
            Adum(iR+1,iC)    = k_c;
            Adum(iR-1,iC)    = k_c;
            b(iN) = -q_wire*Dx^2;
%       
            
        end
        %         figure(2);imagesc(Adum) % These commands can be helpful to debug.
%         set(gca,'YDir','normal')
        
        A(iN,:)=Adum(idx); % convert Adum to vector and place as row in A matrix
end

%b = b(:);

tic
T = A\b; % solves system of equations. Produces temperature as a vector T
toc
Tall=zeros(R,C);
Tall(idx)= T; % reshape T vector to correspond to physical positions.


%% display result

% figure(2);imagesc(A)
% title('geometry')
% colorbar
% set(gca,'YDir','normal')

Tall(Tall==0)=30;
figure(3);imagesc(Tall, [0 80]) % display result ,[5*floor(min(T(:))/5) 5*ceil(max(T(:))/5)]
set(gca,'YDir','normal')
c = colorbar;
c.Label.String = 'Temperature (?C)';
title('Temperature')
colormap hot

            
            
            
            
            

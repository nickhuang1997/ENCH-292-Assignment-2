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
Dx = 0.0005;       %.5 mm cell width and height FINE GRIDe
%Dx = 0.001;         %UNCOMMENT FOR COARSE GRIDE

k_c       = 385;     %thermal conductivity of copper wire   %%%%% engr toolbox
k_i       = 2;     %thermal conductivity of insulation
k_t       = 4;     %thermal conductivity of tile
k_air     = 0.0257; %from engineering toolbox
hc = 19.91          %k*Nu*z (average HTC for forced convection of model)
T_air = 20         %in KELVIN whoohooo

Bi = hc*Dx/k_air


q_wire    = 400/0.000025 %W m-3 
R = int16(h/Dx);
C = w/Dx;       
%%
% initialising the geometry

geom = zeros(R,C);
geom(1:end,1:end) = 1;

idx=find(geom==1);      %finds all the cells with 1's 
                        %(i.e. in the domain. should be all the cells in our case)
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
if Dx == 0.0005
    for iN = 1:length(idx)
        Adum = zeros(R,C);          % Dummy 2D matrix containing weighting at each cell
        Adum(idx(iN)) = 2;          % set current entry to a value of 2
        [iR, iC]=find(Adum==2);     % find locations of the entry identified on above row

    %%%%%%%%%%%%%%%%%%_corners of model_%%%%%%%%%%%%%%%%%%%%

            if iR == 1 && iC == 1   % Bottom LH corner cell
                cas = 1;
                Adum(iR,iC:iC+1) = [-2 1];
                Adum(iR+1,iC)    = 1;

                type(iN) = 1;

            elseif iR == 80 && iC == 1      %top left of model - exposed to air 
                cas = 2;
                Adum(iR,iC:iC+1) = [-(1 + Bi) .5];
                Adum(iR-1,iC) = .5;
                b(iN) = -Bi*T_air

                type(iN) = 2;

            elseif iR == 80 && iC == 66      %top right of model
                cas = 3;
                Adum(iR,iC-1:iC) = [.5 -(1 + Bi)];
                Adum(iR-1,iC) = .5;
                b(iN) = -Bi*T_air;

                type(iN) = 3;

            elseif iR == 1 && iC == 66      %bottom right of model
                cas = 4;
                Adum(iR,iC-1:iC) = [1 -2];
                Adum(iR+1,iC)    = 1;

                type(iN) = 4;
    %         
    %%%%%%%%%%%%%%%%%%_where insulation meets tile_%%%%%%%%%%%%%%%%%%%%
    %                      right and left side                        %



            elseif iR == 40 && iC == 66      %right side insulation meets tile
                cas = 6;
                Adum(iR,iC-1:iC) = [k_t -(2*k_t+k_i)];
                Adum(iR+1,iC)    = k_t;
                Adum(iR-1,iC)    = k_i;

                type(iN) = 5;

            elseif iR == 40 && iC == 1      %right side insulation meets tile
                cas = 6;
                Adum(iR,iC:iC+1) = [-(2*k_t+k_i) k_t];
                Adum(iR+1,iC)    = k_t;
                Adum(iR-1,iC)    = k_i;

                type(iN) = 24;   
                %MAKES BOTTOM HALF MORE YELLOW
            elseif iR == 40 && iC > 1 && iC < 66 %do range for bottom of tile
                cas = 24;
                Adum(iR,iC-1:iC+1) = [k_t -(3*k_t+k_i) k_t];
                Adum(iR+1,iC)    = k_t;
                Adum(iR-1,iC)    = k_i;

                type(iN) = 6;

    %%%%%%%%%%%%%%%%%%_corners of wire _%%%%%%%%%%%%%%%%%%%%

            elseif iR == 10 && iC == 1      %bottom left of wire
                cas = 7;
                Adum(iR,iC:iC+1) = [-(2*k_c+k_i) k_c];
                Adum(iR+1,iC) = k_c;
                Adum(iR-1,iC) = k_i;
                b(iN) = -q_wire*Dx^2;   %done by inspection 

                type(iN) = 7;

            elseif iR == 20 && iC == 1      %top left of wire
                cas = 8;
                Adum(iR,iC:iC+1) = [-(2*k_c+k_i) k_c];
                Adum(iR+1,iC) = k_i;
                Adum(iR-1,iC) = k_c;
                b(iN) = -q_wire*Dx^2;  %done   
                type(iN) = 8;
            elseif iR == 10 && iC == 5      %bottom right corner of the wire
                cas = 9;
                Adum(iR,iC-1:iC+1) = [k_c -(2*k_c+2*k_i) k_i];
                Adum(iR+1,iC)    = k_c;
                Adum(iR-1,iC)    = k_i;
                b(iN) = -q_wire*Dx^2;  %done
                type(iN) = 9;
            elseif iR == 20 && iC == 5      %top right corner of the wire
                cas = 10;
                Adum(iR,iC-1:iC+1) = [k_c -(2*k_c+2*k_i) k_i];
                Adum(iR+1,iC)    = k_i;
                Adum(iR-1,iC)    = k_c;
                b(iN) = -q_wire*Dx^2;  %done

                type(iN) = 10;
    %%%%%%%%%%%%%%%%%%_top and bottom of wire_%%%%%%%%%%%%%%%%%%%%

            elseif iR == 20 && iC < 6 %do range for top of copper
                cas = 11;
                Adum(iR,iC-1:iC+1) = [k_c -(3*k_c+k_i) k_c];
                Adum(iR+1,iC)    = k_i;
                Adum(iR-1,iC)    = k_c;
                b(iN) = -q_wire*Dx^2;  %done

                type(iN) = 11;

            elseif iR == 10 && iC < 6    %do range for bottom of copper
                cas = 12;
                Adum(iR,iC-1:iC+1) = [k_c -(3*k_c+k_i) k_c];
                Adum(iR+1,iC)    = k_c;
                Adum(iR-1,iC)    = k_i;
                b(iN) = -q_wire*Dx^2;  %done    

                type(iN) = 12;
    %%%%%%%%%%%%%%%%%%_top of model_%%%%%%%%%%%%%%%%%%%%

        %RANGE PROBLEMS WITH LINE 170? %
            elseif iR == 80 && iC > 0   %do range for top of t exposed to air
                cas = 13;
                Adum(iR,iC-1:iC+1) = [.5 -(2+Bi) .5];
    %             Adum(iR+1,iC)    = ;
                Adum(iR-1,iC)    = 1;
                b(iN) = -T_air*Bi; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                type(iN) = 13;

    %%%%%%%%%%%%%%%%%%_interior cells_%%%%%%%%%%%%%%%%%%%%

            elseif iR > 10 && iR < 20 && iC < 5 && iC > 1 %range for all interior wire cells 
                cas = 14;
                Adum(iR,iC-1:iC+1) = [1 -(4) 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                b(iN) = -q_wire*Dx^2/k_c;  %done

                type(iN) = 14;

            elseif  iR > 39 && iC > 1 && iC < 66 %range for all interior tile cells 
                cas = 16;
                Adum(iR,iC-1:iC+1) = [1 -(4) 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;

                type(iN) = 16;
    % 
    % %%%%%%%%%%%%%%%%%%_bottom of insulation_%%%%%%%%%%%%%%%%%%%%


        %%%%%%  breaks it ???????????  %%%%%%%%%%%%%%%%%%


            elseif iR < 2 && iC > 1 && iC < 66 %do range for all bottom of the insulation            
                cas = 17;
                Adum(iR,iC-1:iC+1)  = [1 -(3) 1];
                Adum(iR+1,iC)       = 1;

                type(iN) = 17;
    %%%%%%%%%%%%%%%%%%%%%_left walls_%%%%%%%%%%%%%%%%%%%%%%%

            elseif (iR > 1 && iR < 10 && iC < 2) | (iR > 20 && iR < 40 && iC < 2)  %do range for all left walls of insulation
                cas = 18;
                Adum(iR,iC:iC+1) = [-3 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;

                type(iN) = 18;

            elseif iR > 39 && iR < 80 && iC < 2%do range for all left walls of tile
                cas = 19;
                Adum(iR,iC:iC+1) = [-3 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;

                type(iN) = 19;

            elseif iR > 10 && iR < 20 && iC < 2%do range for all left walls of copper
                cas = 20;
                Adum(iR,iC:iC+1) = [-3 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                b(iN) = -q_wire*Dx^2/k_c;
                type(iN) = 20;
    %%%%%%%%%%%%%%%%%%%%%_right walls_%%%%%%%%%%%%%%%%%%%%%%%
    % 
            elseif iR > 1 && iR < 40 && iC > 65 %do range for all right walls of insulation
                cas = 21;
                Adum(iR,iC-1:iC) = [1 -3];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                type(iN) = 21;

            elseif iR > 39 && iR < 80 && iC > 65%do range for all right walls of tile
                cas = 22;
                Adum(iR,iC-1:iC) = [1 -3];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                type(iN) = 22;

            elseif iR < 20 && iR > 9 && iC == 5 %do range for all right walls of copper
                cas = 23;
                Adum(iR,iC-1:iC) = [1 -3];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                b(iN) = -q_wire*Dx^2/k_c;
                type(iN) = 23;


            else%range for all interior insulation cells 

                cas = 15;
                Adum(iR,iC-1:iC+1) = [1 -(4) 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                type(iN) = 15;

            end
            %         figure(2);imagesc(Adum) % These commands can be helpful to debug.
    %         set(gca,'YDir','normal')

            A(iN,:)=Adum(idx); % convert Adum to vector and place as row in A matrix
    end
    
else %goes through the calculations for a COARSE GRID
        for iN = 1:length(idx)
        Adum = zeros(R,C);          % Dummy 2D matrix containing weighting at each cell
        Adum(idx(iN)) = 2;          % set current entry to a value of 2
        [iR, iC]=find(Adum==2);     % find locations of the entry identified on above row

    %%%%%%%%%%%%%%%%%%_corners of model_%%%%%%%%%%%%%%%%%%%%

            if iR == 1 && iC == 1   % Bottom LH corner cell
                cas = 1;
                Adum(iR,iC:iC+1) = [-2 1];
                Adum(iR+1,iC)    = 1;

                type(iN) = 1;

            elseif iR == 40 && iC == 1      %top left of model - exposed to air 
                cas = 2;
                Adum(iR,iC:iC+1) = [-(1 + Bi) .5];
                Adum(iR-1,iC) = .5;
                b(iN) = -Bi*T_air

                type(iN) = 2;

            elseif iR == 40 && iC == 33      %top right of model
                cas = 3;
                Adum(iR,iC-1:iC) = [.5 -(1 + Bi)];
                Adum(iR-1,iC) = .5;
                b(iN) = -Bi*T_air;

                type(iN) = 3;

            elseif iR == 1 && iC == 33      %bottom right of model
                cas = 4;
                Adum(iR,iC-1:iC) = [1 -2];
                Adum(iR+1,iC)    = 1;

                type(iN) = 4;
    %         
    %%%%%%%%%%%%%%%%%%_where insulation meets tile_%%%%%%%%%%%%%%%%%%%%
    %                      right and left side                        %



            elseif iR == 20 && iC == 33      %right side insulation meets tile
                cas = 6;
                Adum(iR,iC-1:iC) = [k_t -(2*k_t+k_i)];
                Adum(iR+1,iC)    = k_t;
                Adum(iR-1,iC)    = k_i;

                type(iN) = 5;

            elseif iR == 20 && iC == 1      %right side insulation meets tile
                cas = 6;
                Adum(iR,iC:iC+1) = [-(2*k_t+k_i) k_t];
                Adum(iR+1,iC)    = k_t;
                Adum(iR-1,iC)    = k_i;

                type(iN) = 24;   
                %MAKES BOTTOM HALF MORE YELLOW
            elseif iR == 20 && iC > 1 && iC < 33 %do range for bottom of tile
                cas = 24;
                Adum(iR,iC-1:iC+1) = [k_t -(3*k_t+k_i) k_t];
                Adum(iR+1,iC)    = k_t;
                Adum(iR-1,iC)    = k_i;

                type(iN) = 6;

    %%%%%%%%%%%%%%%%%%_corners of wire _%%%%%%%%%%%%%%%%%%%%

            elseif iR == 5 && iC == 1      %bottom left of wire
                cas = 7;
                Adum(iR,iC:iC+1) = [-(2*k_c+k_i) k_c];
                Adum(iR+1,iC) = k_c;
                Adum(iR-1,iC) = k_i;
                b(iN) = -q_wire*Dx^2;   %done by inspection 

                type(iN) = 7;

            elseif iR == 10 && iC == 1      %top left of wire
                cas = 8;
                Adum(iR,iC:iC+1) = [-(2*k_c+k_i) k_c];
                Adum(iR+1,iC) = k_i;
                Adum(iR-1,iC) = k_c;
                b(iN) = -q_wire*Dx^2;  %done   
                type(iN) = 8;
            elseif iR == 5 && iC == 3      %bottom right corner of the wire %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cas = 9;
                Adum(iR,iC-1:iC+1) = [k_c -(2*k_c+2*k_i) k_i];
                Adum(iR+1,iC)    = k_c;
                Adum(iR-1,iC)    = k_i;
                b(iN) = -q_wire*Dx^2;  %done
                type(iN) = 9;
            elseif iR == 10 && iC == 3      %top right corner of the wire
                cas = 10;
                Adum(iR,iC-1:iC+1) = [k_c -(2*k_c+2*k_i) k_i];
                Adum(iR+1,iC)    = k_i;
                Adum(iR-1,iC)    = k_c;
                b(iN) = -q_wire*Dx^2;  %done

                type(iN) = 10;
    %%%%%%%%%%%%%%%%%%_top and bottom of wire_%%%%%%%%%%%%%%%%%%%%

            elseif iR == 20 && iC < 6 %do range for top of copper
                cas = 11;
                Adum(iR,iC-1:iC+1) = [k_c -(3*k_c+k_i) k_c];
                Adum(iR+1,iC)    = k_i;
                Adum(iR-1,iC)    = k_c;
                b(iN) = -q_wire*Dx^2;  %done

                type(iN) = 11;

            elseif iR == 10 && iC < 6    %do range for bottom of copper
                cas = 12;
                Adum(iR,iC-1:iC+1) = [k_c -(3*k_c+k_i) k_c];
                Adum(iR+1,iC)    = k_c;
                Adum(iR-1,iC)    = k_i;
                b(iN) = -q_wire*Dx^2;  %done    

                type(iN) = 12;
    %%%%%%%%%%%%%%%%%%_top of model_%%%%%%%%%%%%%%%%%%%%

        %RANGE PROBLEMS WITH LINE 170? %
            elseif iR == 80 && iC > 0   %do range for top of t exposed to air
                cas = 13;
                Adum(iR,iC-1:iC+1) = [.5 -(2+Bi) .5];
    %             Adum(iR+1,iC)    = ;
                Adum(iR-1,iC)    = 1;
                b(iN) = -T_air*Bi; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                type(iN) = 13;

    %%%%%%%%%%%%%%%%%%_interior cells_%%%%%%%%%%%%%%%%%%%%

            elseif iR > 10 && iR < 20 && iC < 5 && iC > 1 %range for all interior wire cells 
                cas = 14;
                Adum(iR,iC-1:iC+1) = [1 -(4) 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                b(iN) = -q_wire*Dx^2/k_c;  %done

                type(iN) = 14;

            elseif  iR > 39 && iC > 1 && iC < 66 %range for all interior tile cells 
                cas = 16;
                Adum(iR,iC-1:iC+1) = [1 -(4) 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;

                type(iN) = 16;
    % 
    % %%%%%%%%%%%%%%%%%%_bottom of insulation_%%%%%%%%%%%%%%%%%%%%


        %%%%%%  breaks it ???????????  %%%%%%%%%%%%%%%%%%


            elseif iR < 2 && iC > 1 && iC < 66 %do range for all bottom of the insulation            
                cas = 17;
                Adum(iR,iC-1:iC+1)  = [1 -(3) 1];
                Adum(iR+1,iC)       = 1;

                type(iN) = 17;
    %%%%%%%%%%%%%%%%%%%%%_left walls_%%%%%%%%%%%%%%%%%%%%%%%

            elseif (iR > 1 && iR < 10 && iC < 2) | (iR > 20 && iR < 40 && iC < 2)  %do range for all left walls of insulation
                cas = 18;
                Adum(iR,iC:iC+1) = [-3 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;

                type(iN) = 18;

            elseif iR > 39 && iR < 80 && iC < 2%do range for all left walls of tile
                cas = 19;
                Adum(iR,iC:iC+1) = [-3 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;

                type(iN) = 19;

            elseif iR > 10 && iR < 20 && iC < 2%do range for all left walls of copper
                cas = 20;
                Adum(iR,iC:iC+1) = [-3 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                b(iN) = -q_wire*Dx^2/k_c;
                type(iN) = 20;
    %%%%%%%%%%%%%%%%%%%%%_right walls_%%%%%%%%%%%%%%%%%%%%%%%
    % 
            elseif iR > 1 && iR < 40 && iC > 65 %do range for all right walls of insulation
                cas = 21;
                Adum(iR,iC-1:iC) = [1 -3];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                type(iN) = 21;

            elseif iR > 39 && iR < 80 && iC > 65%do range for all right walls of tile
                cas = 22;
                Adum(iR,iC-1:iC) = [1 -3];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                type(iN) = 22;

            elseif iR < 20 && iR > 9 && iC == 5 %do range for all right walls of copper
                cas = 23;
                Adum(iR,iC-1:iC) = [1 -3];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                b(iN) = -q_wire*Dx^2/k_c;
                type(iN) = 23;


            else%range for all interior insulation cells 

                cas = 15;
                Adum(iR,iC-1:iC+1) = [1 -(4) 1];
                Adum(iR+1,iC)    = 1;
                Adum(iR-1,iC)    = 1;
                type(iN) = 15;

            end
            %         figure(2);imagesc(Adum) % These commands can be helpful to debug.
    %         set(gca,'YDir','normal')

            A(iN,:)=Adum(idx); % convert Adum to vector and place as row in A matrix
        end
        
end

    
b = b(:);

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
% 
% geo = reshape(type,iR,iC);
% figure(420)
% imagesc(geo)
% set(gca,'YDir','normal')


figure(3);imagesc(Tall, [-5 250]) % display result ,[5*floor(min(T(:))/5) 5*ceil(max(T(:))/5)]
set(gca,'YDir','normal')
c = colorbar;
c.Label.String = 'Temperature (?C)';
title('Temperature')
colormap hot

mirrorT = fliplr(Tall);
newT = [mirrorT Tall];
floor = [newT newT newT newT newT newT newT newT newT newT]; %New array for all 10 wires

figure(2)
imagesc(flipud(floor),[-5 250]) % display result
c = colorbar;
c.Label.String = 'Temperature (?C)';
colormap hot
title('Full Temperature Map of Floor Space')
            
            
            
            
            

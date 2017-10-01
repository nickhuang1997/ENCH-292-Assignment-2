%%Setting up matrices
A = sparse(length(idx),,length(idx));
%dont use 'zeros'
Adum = zeros(R,C); % dummy 2d matric containing weighting at each row and cell index. 
%used to assign weights to correct cells in A matrix when Adum is converted to a vector 
%form it comprises one row of the A matrix

%%

%set number of cells in each dir
R = int16(H/Dx) %Dx = 0.4e-3 H = 0.032
geom = zeros(
geom(1:h/Dx,1:
geom(1:end,1:
idx = find (geom==1)' %finds all cells with values within domain

figure(1);imagesc(geom

T = A/b %solves system of equns
tall = zeros(R,C);
tall(idc) = T

%%
display result
Tall(Tall==0) = 1000
figure(3);imagesc(Tall,[5*floor(min(T(:))/5
est(gca,'Ydir','normal')
c = colorbar;
c.Label.String = 'Temp'


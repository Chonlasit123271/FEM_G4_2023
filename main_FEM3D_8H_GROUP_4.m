close all;clear all; clc;
tic
syms   xi eta phi % r,s,t - local coordinates

nN = 8; %number of node for hexahedral element
E=17000e6;% Modulus N/m^2
nu=0.2;%Poisson

%Calling the exported gmsh file
FEM_3D_COLUMN_12x1
%global node coordinates from GMSH (xyz)
COORD=0.01.*msh.POS;
%Connectivity of nodes in each hexagonal elements
MEMCON = msh.HEXAS;
%Total Number of hexahedral element
NE=length(MEMCON(:,1));

%Support boundary conditions 
% [NODE, fix_x, fix_y, fix_z]
% 1 = fix 
% 0 = not fix
support=[8  1   1   1;
        % 609 1   1   1;
        % 308 1   1   1;
        % 611 1   1   1;
        % 608 1   1   1;
        % 1209 1   1   1;
         % 51  1   1   1;
         4  1   1   1;
         % 48  1   1   1;
        %89  1   1   1;
         % 28  1   1   1;
         6  1   1   1;
         % 49  1   1   1;
         2  1   1   1];

nSupport=size(support,1);

%=======Problem Definition=====%
%Force control          = 'd'
%Displacement control   = 'f'
CONTROL = 'd';

%Load conditions
u_input = readmatrix("u_input.csv");
N_step = length(u_input);
% step = [1];
%Applying displacement [NODE delx dely delz]
%input_disp YES = 1
%input_disp NO = 0

uDisp = zeros(length(COORD)*3);

%start of applying load step

Udisp_result = [];
ffem_result = [];

for step = 1:N_step

%Defining input displacement
%format [ NODE uDisp_x uDisp_y uDisp_z]
%if we want to use cylcic disp from load_apply.mlx

input_disp=[  
            % 7 U_load_list(N_step) 0 0;
            % 52 U_load_list(N_step) 0 0;
            % 3 U_load_list(N_step) 0 0;
            % 38 U_load_list(N_step) 0 0;
            % 90 U_load_list(N_step) 0 0;
            % 18 U_load_list(N_step) 0 0;
            % 5 U_load_list(N_step) 0 0;
            % 50 U_load_list(N_step) 0 0;
             % 1 U_load_list(N_step) 0 0 ];
             7 uDisp(7*3-2)+u_input(step) 0 0;
             5  uDisp(5*3-2)+u_input(step) 0 0];


nLoad_disp=size(input_disp,1);


%Applying displacement [NODE delx dely delz]
%input_disp YES = 1
%input_disp NO = 0
%Applying nodal force [NODE fx fy fz]
nodal_force = [1 3000/4 0 0;
            3 3000/4 0 0;
            5 3000/4 0 0;
            7 3000/4 0 0];
n_nodal_f = length(nodal_force(:,1));


for i = 1:NE
    ELEMENT(i).k = zeros(24,24);
end

disp('Calculaing ELEMENT STIFFNESS MATRIX (ke). Please wait...')
%% CALCULATION OF ELEMENT STIFFNESS MATRIX (ke)
% element #number
for eNo =1:NE
  %storing node coordinates
  xp = zeros(1,length(COORD));
  yp = zeros(1,length(COORD));
  zp = zeros(1,length(COORD));
    for i = 1:length(COORD)
    xp(i)= COORD(i,1); 
    yp(i)= COORD(i,2); 
    zp(i)= COORD(i,3);
    end
% shape functions
syms xi eta phi
N(1)=1/8*((1-xi)*(1-eta)*(1-phi));
N(2)=1/8*((1-xi)*(1-eta)*(1+phi));
N(3)=1/8*((1+xi)*(1-eta)*(1+phi));
N(4)=1/8*((1+xi)*(1-eta)*(1-phi));
N(5)=1/8*((1-xi)*(1+eta)*(1-phi));
N(6)=1/8*((1-xi)*(1+eta)*(1+phi));
N(7)=1/8*((1+xi)*(1+eta)*(1+phi));
N(8)=1/8*((1+xi)*(1+eta)*(1-phi));
nN=length(N);

%% --------------------ELEMENT STiFFNESS(ke)-------------------------------
%Coordinate Mapping
%Mapping
x=0;y=0;z=0;
for i=1:nN
    x=x+xp(MEMCON(eNo,i))*N(i);
    y=y+yp(MEMCON(eNo,i))*N(i);
    z=z+zp(MEMCON(eNo,i))*N(i);
end

%Jacobian
J=[diff(x,xi) diff(y,xi) diff(z,xi);
    diff(x,eta) diff(y,eta) diff(z,eta);
    diff(x,phi) diff(y,phi) diff(z,phi)];

%Differentiating the N matrix to get doMat, dN
%consider using J\[A] instead of inv(J)*[A]
for i=1:nN
    dN(:,i)=J\[diff(N(i),xi);diff(N(i),eta);diff(N(i),phi)];
end
%Arranging the terms to get the B matrix
%B=zeros(3,3*length(N));
for i=1:nN
    j=3*i-2; k=3*i;
B(1:6,j:k)=[dN(1,i), 0      ,0 ;
            0      , dN(2,i),0 ;
            0      , 0      ,  dN(3,i); 
            dN(2,i), dN(1,i), 0;
            0      , dN(3,i) ,dN(2,i);
            dN(3,i) ,   0   ,dN(1,i)];
end

ELEMENT(eNo).B = B;

D = E/((1+nu)*(1-2*nu))* ...
    [1 - nu, nu, nu, 0, 0, 0; 
    nu, 1 - nu, nu, 0, 0, 0;
    nu, nu, 1 - nu, 0, 0, 0;
    0, 0, 0, 1/2 - nu, 0, 0; 
    0, 0, 0, 0, 1/2 - nu, 0;
    0, 0, 0, 0, 0, 1/2 - nu];

%Integrand is
f= B' * D * B * det(J);

%Gauss Points
[xg,wg]=lgwt(2,-1,1);
%Gaussian evaluation of integral [2 Gauss Points]
%ELEMENT(eNo).stiffness=int(int(int(f,xi,-1,1),eta,-1,1),phi,-1,1);
element_stiffness=...
    wg(1)*wg(1)*wg(1)*subs(f,[xi,eta,phi],[xg(1),xg(1),xg(1)]) + ...
    wg(2)*wg(1)*wg(1)*subs(f,[xi,eta,phi],[xg(2),xg(1),xg(1)]) + ...
    wg(2)*wg(1)*wg(2)*subs(f,[xi,eta,phi],[xg(2),xg(1),xg(2)]) + ...
    wg(1)*wg(1)*wg(2)*subs(f,[xi,eta,phi],[xg(1),xg(1),xg(2)]) + ...
    wg(1)*wg(2)*wg(1)*subs(f,[xi,eta,phi],[xg(1),xg(2),xg(1)]) + ...
    wg(2)*wg(2)*wg(1)*subs(f,[xi,eta,phi],[xg(2),xg(2),xg(1)]) + ...
    wg(2)*wg(2)*wg(2)*subs(f,[xi,eta,phi],[xg(2),xg(2),xg(2)]) + ...
    wg(1)*wg(2)*wg(2)*subs(f,[xi,eta,phi],[xg(1),xg(2),xg(2)]);

    ELEMENT(eNo).k = eval(element_stiffness);
end
%End of ke calculation
disp('ELEMENT STIFFNESS(ke) calculation has been finished')
toc

disp('Solving finite element problem. Please wait...')
tic
NDOF = length(COORD)*3;
%% Mapping of local stiffness matrices on to the global stiffness matrices
KG(1:NDOF,1:NDOF)=0;                             %KG be the global stiffness matrix
for i=1:NE
  for lr=1:8                                  %lr is the local row
    for lc=1:8                                 %lc is the local columns
      GRow=MEMCON(i,lr);                         %Global row
      GCol=MEMCON(i,lc);                         %Global column

      KG(GRow*3,GCol*3-2) =KG(GRow*3,GCol*3-2)+ ELEMENT(i).k(lr*3,lc*3-2);        %for W
      KG(GRow*3,GCol*3-1) =KG(GRow*3,GCol*3-1)+ELEMENT(i).k(lr*3,lc*3-1); %for V
      KG(GRow*3,GCol*3)   =KG(GRow*3,GCol*3)+ELEMENT(i).k(lr*3,lc*3); %for U

      KG(GRow*3-1,GCol*3-2) =KG(GRow*3-1,GCol*3-2)+ ELEMENT(i).k(lr*3-1,lc*3-2);        %for W
      KG(GRow*3-1,GCol*3-1) =KG(GRow*3-1,GCol*3-1)+ELEMENT(i).k(lr*3-1,lc*3-1); %for V
      KG(GRow*3-1,GCol*3)   =KG(GRow*3-1,GCol*3)+ELEMENT(i).k(lr*3-1,lc*3); %for U

      KG(GRow*3-2,GCol*3-2) =KG(GRow*3-2,GCol*3-2)+ ELEMENT(i).k(lr*3-2,lc*3-2);        %for W
      KG(GRow*3-2,GCol*3-1) =KG(GRow*3-2,GCol*3-1)+ELEMENT(i).k(lr*3-2,lc*3-1); %for V
      KG(GRow*3-2,GCol*3)   =KG(GRow*3-2,GCol*3)+ELEMENT(i).k(lr*3-2,lc*3); %for U

    end
  end
end

KG_og = KG; %storing original global stiffness matrix

%-------------------------------------------------------------------------%
%%           Apply Support Conditions to the Global Stiffnes Matrix       %
%-------------------------------------------------------------------------%
% This part makes zeros all columns and rows where the supports are except
% the diagonal element of the matrix. Diagonal element set to 1. I choose
% this method because with this way sort of displacement evaulated are not
% changes. And later it is easier to use evaluated values. Only negative
% side of this approach is we have to be careful not to put force where the
% support is fixed.
for i=1:nSupport
    n = support(i,1);
    if (support(i,2) == 1)  %for x direction
        KG(3*n-2,:) = 0;
        KG(:,3*n-2) = 0;
        KG(3*n-2,3*n-2) = 1;
    end
     if (support(i,3) == 1) %for y direction
        KG(3*n-1,:) = 0;
        KG(:,3*n-1) = 0;
        KG(3*n-1,3*n-1) = 1;
    end
    if (support(i,4) == 1) %for z direction
        KG(3*n,:) = 0;
        KG(:,3*n) = 0;
        KG(3*n,3*n) = 1;
    end
end



%-------------------------------------------------------------------------%
%                       Load Vector Computation                           %
%-------------------------------------------------------------------------%
% In this part load vector created. If there is a load vector get the value
% from load matrix. If not not load value set to zero.
nNode=length(COORD);
fvec = zeros(nNode*3,1);

if CONTROL == 'f';
    for i=1:n_nodal_f
    n = nodal_force(i,1);    
    fvec(n*3-2) = nodal_force(i,2);   %for f-x (x-direction);
    fvec(n*3-1) = nodal_force(i,3);   %for f-y (y-direction);
    fvec(n*3) = nodal_force(i,4);     %for f-z (z-direction);    
    end
end


%-------------------------------------------------------------------------%
%                       Displacement Boundary Application                 %          
%-------------------------------------------------------------------------%
% In this part load vector created. If there is a load vector get the value
% from load matrix. If not not load value set to zero.
if CONTROL == 'd';
    for i=1:nLoad_disp
    n = input_disp(i,1); 
    %for u (x-direction);
    if input_disp(i,2) ~= 0
    fvec = fvec - input_disp(i,2).*KG(:,3*n-2);
    fvec(3*n-2) = input_disp(i,2);
    KG(3*n-2,:) = 0;KG(:,3*n-2) = 0;KG(3*n-2,3*n-2) = 1;
    end
    %for v (y-direction);
    if input_disp(i,3) ~= 0
    fvec = fvec - input_disp(i,3).*KG(:,3*n-1);
    fvec(3*n-1) = input_disp(i,3);
    KG(3*n-1,:) = 0;KG(:,3*n-1) = 0;KG(3*n-1,3*n-1) = 1;
    end
    %for w (z-direction);
    if input_disp(i,4) ~= 0
    fvec = fvec - input_disp(i,4).*KG(:,3*n);
    fvec(3*n) = input_disp(i,4);
    KG(3*n,:) = 0;KG(:,3*n) = 0;KG(3*n,3*n) = 1;
    end
    end
end

%-------------------------------------------------------------------------%
%                       Displacement Computation                          %
%-------------------------------------------------------------------------%
%consider using K\[f] instead of inv(K)*[f]

uDisp = KG\fvec;


%Storing change in displacement [delta_x, delta_y, delta_z]
COORDdel = zeros(length(COORD):3);
for i = 1:length(uDisp)/3
COORDdel(i,1) = uDisp(3*i-2);
COORDdel(i,2) = uDisp(3*i-1);
COORDdel(i,3) = uDisp(3*i);
end


ffem = KG_og*uDisp;


Force_apply = zeros(length(COORD):3);
%Storing equivalent applied force [force_x, force_y, force_z]
for i = 1:length(fvec)/3
Force_apply(i,1) = ffem(3*i-2);
Force_apply(i,2) = ffem(3*i-1);
Force_apply(i,3) = ffem(3*i);
end

%strain()=B()*uDisp()
%stress=D()*strain(i)
%Calculation of stresses

%Nodal displacements
for i=1:nNode
    NODE(i).u=[uDisp(3*i-2);uDisp(3*i-1); uDisp(3*i)];
end
%Element displacements

for i=1:NE
    temp=[];
    for j=1:nN
        temp=[temp;
              NODE(MEMCON(i,j)).u];
    end
            ELEMENT(i).u=temp;
end

%Element Stresses
SCP = [-1 -1 1,; 1 -1 1 ; 1 1 1 ; -1 1 1 ; -1 1 1; -1 1 1; -1 1 1; -1 1 1];                % Stress calculation points
for i=1:NE
   uelem = ELEMENT(i).u;
   Belem=ELEMENT(i).B;
   strain=Belem*uelem;
   sigma = D*Belem*uelem;
   for j=1:8
       ELEMENT(i).sigma(:,j) =double(subs(sigma,[xi,eta,phi],SCP(j,:))) ;
       ELEMENT(i).strain(:,j) =double(subs(strain,[xi,eta,phi],SCP(j,:))) ;
   end
end 

%=====================save output as file===================
Udisp_result = [Udisp_result uDisp];
ffem_result = [ffem_result ffem];

disp('Solving calculation has been finished')
toc
end

disp('==============================================')
disp('                     DONE                     ')
disp('==============================================')

%===========================================================
csvwrite('fem_u_result.csv',Udisp_result);
csvwrite('fem_f_result.csv',ffem_result);
%===========================================================




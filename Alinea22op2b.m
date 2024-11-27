% Define system parameters
m_lander = 1803.24;
m_rover = 1025;
m = m_lander; 
J = diag([2508.01 1690.54 3522.33]);
z_cg = -0.75; 
Da = -(0.02*1.05*9.5*5)/m;
I3 = eye(3);
O3 = zeros(3);
O31 = [0;0;0];

alpha1 = -10 * pi / 180;
alpha2 = -10 * pi / 180;
alpha3 = 10 * pi / 180;
alpha4 = 10 * pi / 180;

p1 = [1.5 -1.9 z_cg];
p2 = [-1.5 -1.9 z_cg];
p3 = [-1.5 1.9 z_cg];
p4 = [-1.5 -1.9 z_cg];


A = [O3  I3    O3  O3;  
    O3  Da*I3  O3  O3;
    O3  O3  O3  I3;
    O3  O3  O3  O3];


B = [O31  O31  O31  O31;
    1/m * [0; sin(alpha1); -cos(alpha1)]  1/m * [0; sin(alpha2); -cos(alpha2)]  1/m * [0; sin(alpha3); -cos(alpha3)]  1/m * [0; sin(alpha4); -cos(alpha4)];
    O31  O31  O31  O31;
    J\cross(p1, [0; sin(alpha1); -cos(alpha1)])' ...
    J\cross(p2, [0; sin(alpha2); -cos(alpha2)])' ...
    J\cross(p3, [0; sin(alpha3); -cos(alpha3)])' ...
    J\cross(p4, [0; sin(alpha4); -cos(alpha4)])'];


 

C = [I3, O3, O3, O3;
     O3, I3, O3, O3];
disp(C);


% Check Controllability
ctrb_matrix = ctrb(A, B);
controllability_rank = rank(ctrb_matrix);
if controllability_rank == size(A, 1)
    disp('System is controllable');
else
    disp('System is not controllable');
end

% Check Observability
obsv_matrix = obsv(A, C);
observability_rank = rank(obsv_matrix);
if observability_rank == size(A, 1)
    disp('System is observable');
else
    disp('System is not observable');
end

% Check Stability by computing eigenvalues of A
eigenvalues = eig(A);
disp('Eigenvalues of A:');
disp(eigenvalues);
if all(real(eigenvalues) < 0)
    disp('System is stable');
else
    disp('System is not stable');
end
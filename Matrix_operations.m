%Name: Arvind Purohit

clear

%---- setting up the matrices 

P = input('Submit a matrice in form of [A1 A2; B1 B2;]: ');

%P = [3 2; 4 1;];
%P = [3 -1; -2 4;];
%P = [1 -2; -4 8;];
%P = [2 3; -3 8;];           %Test matrix
%P = [1 3; 4 2;];

%using the hand calculation we can generalize the expression as follows
t = trace(P);               %Components of Quadratic polynomial
d = det(P);

e1 = (t + sqrt(t^2 - 4*d))/2 % first value
e2 = (t - sqrt(t^2 - 4*d))/2 % second value

%----first value ratio

r1 = -P(1,2)/(P(1,1) - e1); %ratio using top equation
r2 = -(P(2,2) - e1)/P(2,1); %ratio using bottom equation

check_e1 = r1 - r2 %Should be zero 

%----second value ratio

r1 = -P(1,2)/(P(1,1) - e2); %ratio using top equation
r2 = -(P(2,2) - e2)/P(2,1); %ratio using bottom equation

check_e2 = r1 - r2 %Should be zero 

%------ eigen vectors

v1 = [-P(1,2); (P(1,1) - e1);];
v2 = [-P(1,2); (P(1,1) - e2);];

if v1(1) < 0 %If first element negative
    v1_fix = -v1/sqrt(v1'*v1)
else
    v1_fix = v1/sqrt(v1'*v1)
end
if v2(1) < 0 %If first element negative
    v2_fix = -v2/sqrt(v2'*v2)
else
    v2_fix = v2/sqrt(v2'*v2)
end

%Checking the Eigenvectors 

check1 = P*v1_fix - e1*v1_fix
check2 = P*v2_fix - e2*v2_fix     %Should be a column vector of zeros

X = [v1 v2]

if abs(e1 - e2) < 1e-6
    disp('The matrix is not diagonalizable.')
else
    D = X^-1*P*X
end


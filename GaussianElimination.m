%#5 Gaussian Elimination
% Vanessa Lama
%Math 403-Exam 2
%The main matrix
A=[-4, 6, -1; 1, -3, 4; 3, -7, -7];
%Matrix B
B=[7; -4; -8];

%Augmented form of matrices A and B
AxB=[A, B];
disp("Original Matrix:"), disp(AxB);

%--------------GAUSSIAN ELIMINATION----------------

n = size(A, 2)
 %(i, i) is always our pivot

for i=1:n
  if((AxB(i, i)) == 0)
      AxB([i,n], :)=AxB([n,i], :); 
  endif
%j lets you go down the column i and eliminate  
  for j=i+1:n
    mult = (AxB(j, i))/(AxB(i, i))
    
    %elimination
    AxB(j, :) = (AxB(j, :)) - (mult * AxB(i, :));
    j = j + 1;
  endfor
  %Changing leading entry to the value of 1
  AxB(i, :) = (AxB(i, :))/(AxB(i, i));
  i = i + 1;
endfor
disp("Matrix after Gaussian Elimination:"), disp(AxB);
     


%-----------------BACKWARD SUBSTITUTION------------------------
x = zeros(n, 1);
for a=n:-1:1
  x(a) = ((AxB(a, n+1))-((AxB(a, a+1:n))*x(a+1:n))) / (AxB(a, a));
endfor
disp("Solution:"), disp(x);
  

  
  


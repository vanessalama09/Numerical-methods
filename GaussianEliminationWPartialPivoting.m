%Main matrices
A =[6, -2, 2, 4; 12, -8, 4, 10; 3, -13, 3, 3; -6, 4, 2, -18];
B =[0; -10; -39; -16];

%Augmented matrix
AxB=[A, B];


disp("Original matrix:"), disp(AxB);
n = size(A, 2);


%-------------------------Gaussian Elimination w/ Partial Pivoting-----------------

%Again (i, i) is the pivot point
for i=1:1:n;  

  %Finging the max value from the pivot column, it was quite time consuming to figure this one out,
  %cause you need max value from column i, row i and further below
  [max_value, index] = max(abs(AxB(i:n, i)))
  %Basic re-arranging to find the index of max in the overall aug. matrix
  real_index= (i-1)+index 

  %Swapping of rows; Partial Pivoting doing it's job
   if(AxB(i,i) != max_value)
     AxB([i, real_index], :) = AxB([real_index, i], :);
   endif
   
  for j=i+1:1:n
    
    
    mult = (AxB(j, i))/(AxB(i, i));
    AxB(j, :) = (AxB(j, :)) - (mult * AxB(i, :));

  endfor
endfor

disp("Matrix after GE w/ PP:"), disp(AxB);


%---------------------------------Backward Substitution-----------------------

x = zeros(n, 1);
for a=n:-1:1
  x(a) = ((AxB(a, n+1))-((AxB(a, a+1:n))*x(a+1:n))) / (AxB(a, a));
endfor
disp(x);

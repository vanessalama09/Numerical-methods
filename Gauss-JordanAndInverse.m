%Main matrix
F = [3, -4, 11, -10; -4, 8 -20, 17; 2, -4, 10, -7; 1, -2, 3, 4]
I = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]
n = size(F, 2)

%Augmented form
FxI=[F, I]

disp("Original Matrix:"), disp(FxI);

%--------------------------The same Gaussian Elimination---------------------
for i=1:n
  if((FxI(i, i)) == 0)
      FxI([i,n], :)=FxI([n,i], :); 
  endif
  
  for j=i+1:n
    mult = (FxI(j, i))/(FxI(i, i));
    FxI(j, :) = (FxI(j, :)) - (mult * FxI(i, :));
    j = j + 1;
  endfor
  FxI(i, :) = (FxI(i, :))/(FxI(i, i));
  i = i + 1;
endfor
disp("Echelon Form:"), disp(FxI);

%----------Gaussian-Jordan Elimination for Reduced Echelon Form---------------
for x=n:-1:1
  for y=x-1:-1:1
    mult1 = (FxI(y, x))/(FxI(x, x));
    FxI(y, :) = (FxI(y, :)) - (mult1 * FxI(x, :));
  endfor
endfor
%No need for Substitution methods!!
disp("Reduced Echelon:");
disp(FxI);

%Creating a new matrix and picking only columns 5, 6, 7, and 8 of FxI
Finv = FxI(:, [5 6 7 8]);

disp("Inverse of F:"), disp(Finv);

%-------Multiplication of matrix F and Finv
disp("To show that multiplying F with Finv really gives the Identity matrix:"), disp(F*Finv);
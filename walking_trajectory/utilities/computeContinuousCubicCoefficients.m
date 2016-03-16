function coefficients = computeContinuousCubicCoefficients(xKnots, durations, xDotInitial, xDotFinal)
   
   if nargin == 2
     xDotInitial = 0;
     xDotFinal = 0;
   end

   numberOfKnots = length(xKnots)

   rowIndex = 0;
   colIndex = 0;
   if numberOfKnots == 2
   	  c(1, 1) = 1;
   	  c(2, 2) = 1;
   	  c(3, 1:4) = [1 durations(1) durations(1)^2 durations(1)^3];
   	  c(4, 1:4) = [0 1 2*durations(1) 3 * durations(1)^2];

   	  b(1) = xKnots(1);
   	  b(2) = xDotInitial;
   	  b(3) = xKnots(2);
   	  b(4) = xDotFinal;
   else
     for i = 1:numberOfKnots-1
     	  if i == 1
     	  	c(1, 1) = 1;
    	    c(2, 2) = 1;
    	    c(3, 1:4) = [1 durations(1) durations(1)^2 durations(1)^3];
   	        c(4, 2:6) = [1 2*durations(1) 3*durations(1)^2 0 -1];
   	        c(5, 3:7) = [2 6*durations(1) 0 0 -1];

   	        b(1) = xKnots(i);
   	        b(2) = xDotInitial;
   	        b(3) = xKnots(i+1);

    	    rowIndex = 5;
    	    colIndex = 4;
    	  elseif i == numberOfKnots-1
    	  	c(rowIndex+1, colIndex+1) = 1;
    	    c(rowIndex+2, colIndex+1:colIndex+4) = [1 durations(i) durations(i)^2 durations(1)^3];
    	    c(rowIndex+3, colIndex+2:colIndex+4) = [1 2*durations(i) 3*durations(i)^2];

    	    b(rowIndex+1) = xKnots(i);
    	    b(rowIndex+2) = xKnots(i+1);
    	    b(rowIndex+3) = xDotFinal;
    	  else
    	  	c(rowIndex+1, colIndex+1) = 1;
    	    c(rowIndex+2, colIndex+1:colIndex+4) = [1 durations(i) durations(i)^2 durations(i)^3];
    	    c(rowIndex+3, colIndex+2:colIndex+6) = [1 2*durations(i) 3*durations(i)^2 0 -1];
    	    c(rowIndex+4, colIndex+3:colIndex+7) = [2 6*durations(i) 0 0 -1];

    	    b(rowIndex+1) = xKnots(i);
    	    b(rowIndex+2) = xKnots(i+1);

    	    rowIndex = rowIndex + 4;
    	    colIndex = colIndex + 4;
    	  end
      end
   end

   coefficients_temp = c^-1 * b';

   rowIndex = 0;
   for i = 1:numberOfKnots-1
   	  coefficients{i} = [coefficients_temp(rowIndex+1:rowIndex+4)];
   	  rowIndex = rowIndex + 4;
   end
end
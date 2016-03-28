function [coefficients, coefficientsDot, coefficientsDotDot] = ...
    computeCubicProjection(durations, numberOfKnots)
    
   if nargin == 1
       numberOfKnots = length(durations) + 1;
   end
   
   if (length(durations) ~= numberOfKnots-1)
       disp('There is an error in the dimensionality of the problem')
   end


      B = zeros(4*(numberOfKnots-1), numberOfKnots);
     
     rowIndex = 0;
     colIndex = 0;

     if numberOfKnots == 2
   	   c(1, 1) = 1;
   	   c(2, 2) = 1;
   	   c(3, 1:4) = [1 durations(1) durations(1)^2 durations(1)^3];
   	   c(4, 1:4) = [0 1 2*durations(1) 3 * durations(1)^2];
       
       B(1,1) = 1;
   	   B(3,2) = 1;
     else
        for i = 1:numberOfKnots-1
     	    if i == 1
     	  	  c(1, 1) = 1;
    	      c(2, 2) = 1;
    	      c(3, 1:4) = [1 durations(1) durations(1)^2 durations(1)^3];
   	          c(4, 2:6) = [1 2*durations(1) 3*durations(1)^2 0 -1];
   	          c(5, 3:7) = [2 6*durations(1) 0 0 -1];
         
              B(1,i) = 1;
              B(3,i+1) = 1;
              
    	      rowIndex = 5;
    	      colIndex = 4;
    	    elseif i == numberOfKnots-1
    	  	  c(rowIndex+1, colIndex+1) = 1;
    	      c(rowIndex+2, colIndex+1:colIndex+4) = [1 durations(i) durations(i)^2 durations(1)^3];
    	      c(rowIndex+3, colIndex+2:colIndex+4) = [1 2*durations(i) 3*durations(i)^2];
              
              B(rowIndex+1,i) = 1;
              B(rowIndex+2,i+1) = 1;
    	    else
    	      c(rowIndex+1, colIndex+1) = 1;
    	      c(rowIndex+2, colIndex+1:colIndex+4) = [1 durations(i) durations(i)^2 durations(i)^3];
    	      c(rowIndex+3, colIndex+2:colIndex+6) = [1 2*durations(i) 3*durations(i)^2 0 -1];
    	      c(rowIndex+4, colIndex+3:colIndex+7) = [2 6*durations(i) 0 0 -1];
              
              B(rowIndex+1,i) = 1;
              B(rowIndex+2,i+1) = 1;
              
    	      rowIndex = rowIndex + 4;
    	      colIndex = colIndex + 4;
    	    end
        end
     end

     coefficients_temp = c^-1 * B;
     
     rowIndex = 0;
     for i = 1:numberOfKnots-1
   	    coefficients{i} = [coefficients_temp(rowIndex+1:rowIndex+4,:)];
        
        coefficientsDot{i}(1,:) = coefficients{i}(2,:);
        coefficientsDot{i}(2,:) = 2 * coefficients{i}(3,:);
        coefficientsDot{i}(3,:) = 3 * coefficients{i}(4,:);
        
        coefficientsDotDot{i}(1,:) = coefficientsDot{i}(2,:);
        coefficientsDotDot{i}(2,:) = 2 * coefficientsDot{i}(3,:);
        
   	    rowIndex = rowIndex + 4;
     end
end
%EVALUATEFOURIERSERIES Returns a vector of phi values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function takes a series of coefficients and a vector of times.  It
%then returnes a vector of phi values that correspond to the values in the
%time array it was passed.  The format of the coefficients must be
%initial,a1,b1,a2,b2...an,bn where a corresponds to sin terms and b to cos.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by PlotAndWriteToFile
function PhiPoints = EvaluateFourierSeries(direction,c,Tarray)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function calculates phi at a given time value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by PhiPoints
    function phi = Phi(c,t)
        
        %The fourier series starts with one coefficient at order 0 and then
        %gains two coeffiecients with each increase in order.
        order = (length(c)-1)/2;
        %The first element is divided by two and does not follow the normal
        %pattern
        phi=c(1)/2;
        
        %Gor every order afterwards there is a sin and cos term
        for j=1:order 
            %The following provides the ability to run the two trips with
            %different fourier frequencies.
            if direction == Direction.FORWARD
                phi=phi+c(2*j)*cos(j*t/2)+c(2*j+1)*sin(j*t/2);
            elseif direction == Direction.RETURN
                phi=phi+c(2*j)*cos(j*t/2)+c(2*j+1)*sin(j*t/2);
            end
        end
    end
%fills the PhiPoints vector with values corresponding to the Time vector
for i = 1:length(Tarray)
    PhiPoints(i) = Phi(c,Tarray(i));
end
end

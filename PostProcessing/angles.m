function [gamma,tau] = angles(time,tFinal, gammaCoeff,tauCoeff)
    time = time/tFinal;

    gamma = gammaCoeff(1);
    for m=1:(length(gammaCoeff)-1)/2
        gamma=gamma+gammaCoeff(2*m)*cos(2*pi*m*time)+gammaCoeff(2*m+1)*sin(2*pi*m*time);
    end
    tau = tauCoeff(1);
    for m=1:(length(tauCoeff)-1)/2
        tau=tau+tauCoeff(2*m)*cos(2*pi*m*time)+tauCoeff(2*m+1)*sin(2*pi*m*time);
    end
end


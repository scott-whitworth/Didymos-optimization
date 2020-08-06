%STORELOWESTVALUES Store values for lowest converging time guess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is passed all variables that need to be changed if a new
%convergence is found for a lower time guess.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Called by MinimizeTripTime
function [lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,lowestYdata,lowestTdata] =...
    StoreLowestValues(converge,lastConverge,cConverge,lastEscape,lastfval,finalMass,escapeVelocity,escapeEarthAngle,earthConditions,asteroidConditions,departBound,approachBound,coastFraction,lowestConverge,lowestcConverge,lowestEscape,lowestfval,lowestfinalMass,lowestescapeVelocity,lowestescapeEarthAngle,lowestearthConditions,lowestasteroidConditions,lowestdepartBound,lowestapproachBound,lowestcoastFraction,flagStore,lowestYdata,lowestTdata,lastYdata,lastTdata)

if (converge && lowestConverge > lastConverge) && flagStore
    lowestConverge = lastConverge;
    lowestcConverge = cConverge;
    lowestEscape = lastEscape;
    lowestfval = lastfval;
    lowestfinalMass = finalMass;
    lowestescapeVelocity = escapeVelocity;
    lowestescapeEarthAngle = escapeEarthAngle;
    lowestearthConditions = earthConditions;
    lowestasteroidConditions = asteroidConditions;
    lowestdepartBound = departBound;
    lowestapproachBound = approachBound;
    lowestcoastFraction = coastFraction;
    lowestTdata = lastTdata;
    lowestYdata = lastYdata;
end
end


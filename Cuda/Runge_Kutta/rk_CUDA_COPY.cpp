


void rk4sys_C_CLONE(Individual *individuals, double *timeInitial, double *startStepSize, double *absTolInput,
                                                  int n, thruster<double> thrust, const cudaConstants* cConstant){
    int threadId = 0;
    if (threadId < n) {
        rkParameters<double> threadRKParameters = individuals[threadId].startParams; // get the parameters for this thread

        elements<double> curPos = threadRKParameters.y0; // start with the initial conditions of the spacecraft

        // storing copies of the input values
        double stepSize = *startStepSize;
        double absTol = *absTolInput;
        double curTime = *timeInitial;
        double startTime = *timeInitial;
        double curAccel = 0;

        elements<double> k1, k2, k3, k4, k5, k6, k7; // k variables for Runge-Kutta calculation of y based off the spacecraft's final state

        double massFuelSpent = 0; // mass of total fuel expended (kg) starts at 0

        bool coast; // to hold the result from calc_coast()

        elements<double> error; // holds output of previous value from rkCalc

        double deltaT = stepSize;

        double n = 0.0;

        while (curTime < threadRKParameters.tripTime) {
            deltaT = stepSize;

            coast = calc_coast(threadRKParameters.coeff, curTime, threadRKParameters.tripTime, thrust);

            curAccel = calc_accel( curPos.r,    curPos.z, thrust, massFuelSpent, deltaT, coast, cConstant->wet_mass, cConstant);

            // calculate k values and get new value of y
            rkCalc(curTime, threadRKParameters.tripTime, stepSize, curPos, threadRKParameters.coeff, curAccel, error, k1, k2, k3, k4, k5, k6, k7, thrust); 

            curTime += stepSize; // update the current time in the simulation
            
            stepSize *= calc_scalingFactor(curPos-error,error,absTol); // Alter the step size for the next iteration

            // The step size cannot exceed the total time divided by 2 and cannot be smaller than the total time divided by 1000
            if (stepSize > (threadRKParameters.tripTime - startTime) / 100) {
                stepSize = (threadRKParameters.tripTime - startTime) / 100;
            }
            else if (stepSize < (threadRKParameters.tripTime - startTime) / 1000) {
                stepSize = (threadRKParameters.tripTime - startTime) / 1000;
            }
            
            if ( (curTime + stepSize) > threadRKParameters.tripTime) {
                stepSize = (threadRKParameters.tripTime - curTime); // shorten the last step to end exactly at time final
            }

            // if the spacecraft is within 0.5 au of the sun, the radial position of the spacecraft artificially increases to 1000, to force that path to not be used in the optimization.
            if ( sqrt(pow(curPos.r,2) + pow(curPos.z,2)) < 0.5) {
                curPos.r = 1000;

                // output to this thread's index
                individuals[threadId].finalPos = curPos;
                individuals[threadId].posDiff = 1.0e10;
                individuals[threadId].velDiff =  0.0;

                return;
            }
            n += 1.0;
        }

         // output to this thread's index
        individuals[threadId].finalPos = curPos;

        //individuals[threadId].getPosDiff(cConstant);
        individuals[threadId].posDiff = n;
        individuals[threadId].getVelDiff(cConstant);
        //individuals[threadId].velDiff = *absTolInput;
        //individuals[threadId].getCost(cConstant);

        return;
    }
    return;
}
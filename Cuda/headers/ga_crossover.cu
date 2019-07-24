// Collection of functions to perform crossover operations on rkParameters
// A crossover mask is an array of elements equal to the number of elements in the rkParameters list
// The mask defines what elements come from partent n and parent m
// [ {elements}     {coefficients}          {other parameters} ]
// [   0-5          6-14,15-19,20-24, 25,    26, 27                                             ]

#include "rkParameters.h"
#include "ga_crossover.h"
#include <iostream>
#include <random>
#include <chrono>

using namespace std;

// Creates a random bifurcation mask
// Randomly picks one index to be the start of the '2's from mask
// input: rng - a constructed mt19937_64 random number generator
// in/out: mask - all data will be overwritten
//              - Based on random index, first selection will be 1's, last selection will be 2's
//              - ex: [1, 1, 1, 1, 2, 2]
void crossOver_randHalf(int mask[], mt19937_64 & rng){
    int crossIndex = rng() % (RK_SIZE-1);
    cout << "Random Index: " << crossIndex << endl;
    for(int i = 0; i < RK_SIZE; i++){
        if(i > crossIndex){
            mask[i] = 2;
        } else {
            mask[i] = 1;
        }        
    }
    return;
}

// Creates a random mask
// Each element in a mask is randomly set to either 1 or 2
// in/out : All data overwritten, set randomly
// input: rng - a constructed mt19937_64 random number generator
void crossOver_wholeRandom(int mask[], mt19937_64 & rng){
    for(int i = 0; i < RK_SIZE; i++ ){
        if(rng()%2){ //Coin flip, either 1/0
            mask[i] = 2;
        } else {
            mask[i] = 1;
        }
    }
    return;
}

//Create a mask the flips just the gamma coefficients
// input mask: over writes all data with [1 ... 1, 2, 2, ... 2, 2, 1 ... 1]
// 2's correspond to Gamma coefficients
void crossOver_gammaPos(int mask[]){
    for(int i = 0; i < RK_SIZE; i++){
        if((i >= 6) && (i <= 14 ) ){
            mask[i] = 2;
        } else {
            mask[i] = 1;
        }
    }
    return;
}

//Create a mask the flips just the tau coefficients
// input mask: over writes all data with [1 ... 1, 2, 2, ... 2, 2, 1 ... 1]
// 2's correspond to tau coefficients
void crossOver_tauPos(int mask[]){
    for(int i = 0; i < RK_SIZE; i++){
        if((i >= 15) && (i <= 19 ) ){
            mask[i] = 2;
        } else {
            mask[i] = 1;
        }
    }
    return;
}

//Utility to flip the polarity of a mask
// input:  mask is an array of size RK_SIZE, input with either 1 or 2 as a mask
// output: each 1 in mask will be reassigned to be a 2, each 2 will be reassigned 1
void flipMask(int mask[]){
    for(int i = 0; i < RK_SIZE; i++){
        if(mask[i] == 1){
            mask[i] = 2;
        } else {
            mask[i] = 1;
        }
    }
    return;
}

//Copy contents of maskIn into maskOut
void copyMask(int maskIn[], int maskOut[]){
    for(int i = 0; i < RK_SIZE; i++){
        maskOut[i] = maskIn[i];
    }
}

void printMask(int mask[]){
    cout << "[";
    for(int i = 0; i < RK_SIZE; i++){
        cout << mask[i];
        if(i < RK_SIZE-1){
            cout <<", ";
        }
    }
    cout << "]";
}

rkParameters<double> generateNewIndividual(const rkParameters<double> & p1, const rkParameters<double> & p2, const int mask[]){
    rkParameters<double> newInd = p1;

    // itterate over set of make values
    for(int i = 0; i < RK_SIZE; i++){
        if(mask[i] == 2){
            //Check elemnents            
            switch(i){
                case 0: //element[0] = r
                    newInd.y0.r = p2.y0.r;
                    break;
                case 1: //element[11 = theta
                    newInd.y0.theta = p2.y0.theta;
                    break;
                case 2: //element[2] = z
                    newInd.y0.z = p2.y0.z;
                    break;
                case 3: //element[3] = v_r
                    newInd.y0.vr = p2.y0.vr;
                    break;                
                case 4: //element[4] = v_theta
                    newInd.y0.vtheta = p2.y0.vtheta;
                    break;                
                case 5: //element[5] = v_z
                    newInd.y0.vz = p2.y0.vz;
                    break;                
            }
            //check coeff
            if( (i > 5) && (i <26) ){
                if(i < 15) {//Gamma (6-14)
                    newInd.coeff.gamma[i-6] = p2.coeff.gamma[i-6];
                } else if(i < 20) {//Tau (15-19)
                    newInd.coeff.tau[i-15] = p2.coeff.tau[i-15];
                } else if(i < 25) {//Coasting (20-24)
                    newInd.coeff.coast[i-20] = p2.coeff.coast[i-20];
                } else if(i < 26) {//Coasting Threshold (25)
                    newInd.coeff.coastThreshold = p2.coeff.coastThreshold;                                   
                } else {
                    //ERROR
                }
            }
            if(i == 26){ //Wetmass
                newInd.wetMass = p2.wetMass;
            }
            if(i == 27){ //Time final
                newInd.timeFinal = p2.timeFinal;
            }
        }
    }

    return newInd;    
}

//Creates a new rkParameter based on the average between p1 and p2
// input: p1 and p2 are valid rkParameters
// output: average of the two
rkParameters<double> generateNewIndividual_avg(const rkParameters<double> & p1, const rkParameters<double> & p2){
    rkParameters<double> newInd;

    newInd.y0.r         = (p1.y0.r/2.0)      + (p2.y0.r/2.0);
    newInd.y0.theta     = (p1.y0.theta/2.0)  + (p2.y0.theta/2.0);
    newInd.y0.z         = (p1.y0.z/2.0)      + (p2.y0.z/2.0);
    newInd.y0.vr        = (p1.y0.vr/2.0)     + (p2.y0.vr/2.0);
    newInd.y0.vtheta    = (p1.y0.vtheta/2.0) + (p2.y0.vtheta/2.0);
    newInd.y0.vz        = (p1.y0.vz/2.0)     + (p2.y0.vz/2.0);

    for(int i = 0; i < p1.coeff.gammaSize; i++){
        newInd.coeff.gamma[i] = (p1.coeff.gamma[i]/2.0) + (p2.coeff.gamma[i]/2.0);
    }
    for(int i = 0; i < p1.coeff.tauSize; i++){
        newInd.coeff.tau[i] = (p1.coeff.tau[i]/2.0) + (p2.coeff.tau[i]/2.0);
    }
    for(int i = 0; i < p1.coeff.coastSize; i++){
        newInd.coeff.coast[i] = (p1.coeff.coast[i]/2.0) + (p2.coeff.coast[i]/2.0);
    }
    newInd.coeff.coastThreshold = (p1.coeff.coastThreshold/2.0) + (p2.coeff.coastThreshold/2.0);
    newInd.wetMass = (p1.wetMass/2.0) + (p2.wetMass/2.0);
    newInd.timeFinal = (p1.timeFinal/2.0) + (p2.timeFinal/2.0);

    return newInd;    
}

void crossover(Individual *survivors, Individual *pool, int selectionSize, int poolSize){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

    int mask[RK_SIZE];
    
    //get masks for every crossover
    for(int i = 0; i < (selectionSize / 4); i++){
        for(int j = 0; j < 4; j++){
            crossOver_wholeRandom(mask, rng);
            pool[poolSize - 1 - (4 * i) - j] = Individual();
            pool[poolSize - 1 - (4 * i) - j].startParams = generateNewIndividual(survivors[2*i].startParams, survivors[(2*i)+1].startParams, mask);
        }
    }
}

//Unit Test for ga_crossover
/*
int main(){
    mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
    
    // input parameters for rk4Simple which are the same for each thread
    coefficients<double> testcoeff;
    for(int i = 0; i < testcoeff.gammaSize; i++){
        testcoeff.gamma[i] = 1;
    }
    for(int i = 0; i < testcoeff.tauSize; i++){
        testcoeff.tau[i] = 1;
    }
    for(int i = 0; i < testcoeff.coastSize; i++){
        testcoeff.coast[i] = 1;
    }
    testcoeff.coastThreshold = 1;
    elements<double> spaceTest(1, 1, 1, 1, 1, 1);
    rkParameters<double> test1(1, 1, spaceTest, testcoeff); 

    // input parameters for rk4Simple which are the same for each thread
    coefficients<double> testcoeff2;
    for(int i = 0; i < testcoeff2.gammaSize; i++){
        testcoeff2.gamma[i] = 2;
    }
    for(int i = 0; i < testcoeff2.tauSize; i++){
        testcoeff2.tau[i] = 2;
    }
    for(int i = 0; i < testcoeff2.coastSize; i++){
        testcoeff2.coast[i] = 2;
    }
    testcoeff2.coastThreshold = 2;
    elements<double> spaceTest2(2, 2, 2, 2, 2, 2);
    rkParameters<double> test2(2, 2, spaceTest2, testcoeff2); 

    cout << "****** Test1 ******\n" << test1 << endl << endl;
    cout << "****** Test2 ******" << test2 << endl;

    //Set Up Mask
    int test_mask[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int test_maskOther[RK_SIZE];
    cout << "               |  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7" << endl;
    cout << "Starting Mask  : ";
    printMask(test_mask);
    cout << endl;

    //crossOver_randHalf(test_mask,rng);
    crossOver_wholeRandom(test_mask,rng);
    copyMask(test_mask,test_maskOther);
    flipMask(test_maskOther);

    cout << "First randHalf : ";
    printMask(test_mask);
    cout << endl;

    cout << "Second randHalf: ";
    printMask(test_maskOther);
    cout << endl;

    //Generating Offspring:
    rkParameters<double> output_1 = generateNewIndividual(test1,test2,test_mask);
    rkParameters<double> output_2 = generateNewIndividual(test1,test2,test_maskOther);
    rkParameters<double> output_3 = generateNewIndividual_avg(test1,test2);

    cout << "****** output_1 ******\n" << output_1 << endl << endl;
    cout << "****** output_2 ******\n" << output_2 << endl << endl;
    cout << "****** output_3 ******\n" << output_3 << endl << endl;



}
*/
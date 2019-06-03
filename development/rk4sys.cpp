#include <vector>

/* fourth-order RUnge-Kutta for a system of ODEs
-integrates a system of ODEs with fourth-order RK method

Input:
tSpan: [ti, tf]; initial and final times with output generated at interval of h,
or [t0 t1 ... tf]; specific times where solution output
y0: initial values of dependent variables
stepSize: the difference between times
output:
tp: independent variables- times
returns: solutions of dependent variables
*/ 
std::vector<std::array<double, 6>> rk4sys(std::vector<double> tspan, double y0[6], double stepSize, std::vector<double> *tp){
    // fill a vector with all the times
    std::vector<double> times;

    double ti = tspan.front();
    double tf = tspan.back();
    double numTimes = tspan.size();
    
    if(numTimes == 2){
        double nextTime = ti;

        while(nextTime < tf){ 
            times.push_back(nextTime);
            nextTime += stepSize; //increment the time by the step amount
        }

        //ensure that the final time is included in case tf is not a multiple of stepSize
        if(times.back() < tf){
            times.push_back(tf);
        }    

        numTimes = times.size();   
    }
    else{
        times = tspan;
    }

    //initialize vectors for the solution
    std::vector<double> t; //store the values to be output in tp
    double tt = ti;

    std::vector<std::array<double, 6>> y;
    //y.push_back(y0);
    for(int i = 0; i < 6; i++){
        y[0][i] = y0[i];
    }
    
    int np = 0;
    t.push_back(tt);
    
    std::vector<std::array<double, 6>> yp;
    yp.push_back(y[0]);
    int j = 0;

    

    while(true){
        double tend = times[np + 1];
        double hh = times[np + 1] - times[np];
        
        if(hh > stepSize)
            hh = stepSize;

            while(true){
                if((tt + hh) > tend)
                    hh = tend - tt;

                //MATH
                double k1, k2, k3, k4;
                //k1 =





                double phi = (k1 + 2 * (k2 + k3) + k4) / 6;
                for(int i = 0; i < 6; i++){
                    y[j+1][i] = y[j][i] + phi * hh;
                }

                //counter
                tt = tt + hh;
                j++;

                if (tt >= tend)
                    break;
            }

            np++;
            t[np] = tt;
            yp[np] = y[j];

            if(tt >= tf)
                break;
    }
    *tp = t;
    return yp;
}

int main(){

}
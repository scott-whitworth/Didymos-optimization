#include <iostream>
#include "d_trip.h"
#include "optimizationCoeff.h"
#include "dataSet.h"


using namespace std;

__global__ void kernel_rk(DataSet * dataSets, D_Trip * output_pos_vel, int n, double finalTime, double stepSize, int loop);

int main(){

    int n = 300;
    double finalTime = 3600*5;
    double stepSize = 5;
    //int iter =850;
    //int iter =900;

    //Host Pointers and Allocation
    DataSet * h_set = new DataSet[n];
    D_Trip * h_out_pv = new D_Trip[n];

    //Fill host Data
    for(int i = 0; i < n; i++){
        h_set[i].pos_vel.v_z = i;
        h_set[i].coeff.gamma[2] = i*2;
    }

    //Event Creation
    cudaEvent_t start_event, stop_event;
    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);
    float k_time = 0;


    //Device Pointers
    DataSet * d_set;
    D_Trip * d_out_pv;

    //Allocate Device Memory
    cudaMalloc((void **) &d_set, n * sizeof(DataSet));
    cudaMalloc((void **) &d_out_pv, n * sizeof(D_Trip));


    //Copy Host to Device
    cudaMemcpy(d_set, h_set, n * sizeof(DataSet), cudaMemcpyHostToDevice);

    cout << "Done Setting up, allocating and transfering data" << endl;

    for(int iter = 100; iter < 1500; iter += 50){

        //Call kernel
        cudaEventRecord(start_event);
        kernel_rk<<<(n+31)/32,32>>>(d_set,d_out_pv,n,finalTime,stepSize,iter);
        cudaEventRecord(stop_event);

        cudaEventSynchronize(stop_event);
        cudaEventElapsedTime(&k_time, start_event, stop_event);


        cout << "Done calling kernel" << endl;
        cout << "#### Loops: \t" << iter << "\tKernel Time: " << k_time << " ms ####" << endl;
        //cout << "#### N: \t" << n << "\tKernel Time: " << k_time << " ms ####" << endl;
        //}
        //Copy off Device
        cudaMemcpy(h_out_pv, d_out_pv, n * sizeof(D_Trip) , cudaMemcpyDeviceToHost);

        cout << "Starting the analyze the data" << endl;
        //Analyze Data
        double testMult = h_out_pv[0].v_r;
        for(int i = 1; i < n; i++){
            if(h_out_pv[i].v_z != h_out_pv[i].v_theta){
                cout << "Error: " << h_out_pv[i] << endl;
                return 0;
            }
            if(testMult != h_out_pv[i].v_r){
                cout << "Error: " << h_out_pv[i] << endl;
                return 0;
            }
        }
        cout << "Last Element: " << h_out_pv[99].v_z << endl;
    }

    //cout << "Last Element: " << h_out_pv[n] << endl;

    cudaFree(d_set);
    cudaFree(d_out_pv);
    
    delete [] h_set;
    delete [] h_out_pv;

    cout << "Done" << endl;
}

__global__ void kernel_rk(DataSet * dataSets, D_Trip * output_pos_vel, int n, double finalTime, double stepSize, int loop){
    int tID = threadIdx.x + (blockIdx.x * blockDim.x);
    if(tID >= n){
        return;
    } else {
        double curTime = 0;
        double mult = 1;
        int loopNum = loop;

        while(curTime < finalTime){
            //Do stuff
            for(int i = 0; i < loopNum; i++){
                mult += pow(curTime,2);
            }

            curTime += stepSize;
        }

        output_pos_vel[tID].v_r = mult;
        output_pos_vel[tID].v_z = dataSets[tID].pos_vel.v_z;
        output_pos_vel[tID].v_theta = tID;
        return;
    }

}
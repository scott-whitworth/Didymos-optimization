

__global__ void kernel_rk(dataSet * dataSets, D_Trip * output_pos_vel, int n, double finalTime, double stepSize){
    int tID = threadIDx.x + (blockIDx.x * blockDim.x);
    if(tID >= n){
        return;
    } else {
        output_pos_vel[tID]->v_theta = tID;
        return;
    }

}
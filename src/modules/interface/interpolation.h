int get_location_dim(int nx, float x[], float xP);
float LI_4D(int nx,
             int ny,
             int nz,
             int nw,
             float x[],
             float y[],
             float z[],
             float w[],
             int f[],
             float xP,
             float yP,
             float zP,
             float wP);
int lower_bound(float arr[], int N, float X);


int get_location_dim(int nx, float x[], float xP) {
    if(xP <= x[0])
        return 0;
    else if(xP >= x[nx-1])
        return nx-2;
    else
        return lower_bound(x, nx, xP);
}

float LI_4D(int nx,
             int ny,
             int nz,
             int nw,
             float x[],
             float y[],
             float z[],
             float w[],
             int f[],
             float xP,
             float yP,
             float zP,
             float wP){

    int i = get_location_dim(nx,x,xP);
    int j = get_location_dim(ny,y,yP);
    int k = get_location_dim(nz,z,zP);
    int l = get_location_dim(nw,w,wP);

    int ip = i+1;
    int jp = j+1;
    int kp = k+1;
    int lp = l+1;

    float xWta = (xP - x[i])/(x[i+1]-x[i]);
    float yWta = (yP - y[j])/(y[j+1]-y[j]);
    float zWta = (zP - z[k])/(z[k+1]-z[k]);
    float wWta = (wP - w[l])/(w[l+1]-w[l]);

    float xWtb = 1.0f-xWta;
    float yWtb = 1.0f-yWta;
    float zWtb = 1.0f-zWta;
    float wWtb = 1.0f-wWta;

    int nwzy = nw*nz*ny;
    int nwz  = nw*nz;

    return ((((float)f[i *nwzy + j *nwz + k *nw + l ]/10000.0f*xWtb +                                // return f[i ][j ][k ][l ]*xWtb*yWtb*zWtb*wWtb +
              (float)f[ip*nwzy + j *nwz + k *nw + l ]/10000.0f*xWta) * yWtb +                        //        f[ip][j ][k ][l ]*xWta*yWtb*zWtb*wWtb +
             ((float)f[i *nwzy + jp*nwz + k *nw + l ]/10000.0f*xWtb +                                //        f[i ][jp][k ][l ]*xWtb*yWta*zWtb*wWtb +
              (float)f[ip*nwzy + jp*nwz + k *nw + l ]/10000.0f*xWta) * yWta) * zWtb +                //        f[ip][jp][k ][l ]*xWta*yWta*zWtb*wWtb +
            (((float)f[i *nwzy + j *nwz + kp*nw + l ]/10000.0f*xWtb +                                //        f[i ][j ][kp][l ]*xWtb*yWtb*zWta*wWtb +
              (float)f[ip*nwzy + j *nwz + kp*nw + l ]/10000.0f*xWta) * yWtb +                        //        f[ip][j ][kp][l ]*xWta*yWtb*zWta*wWtb +    
             ((float)f[i *nwzy + jp*nwz + kp*nw + l ]/10000.0f*xWtb +                                //        f[i ][jp][kp][l ]*xWtb*yWta*zWta*wWtb +
              (float)f[ip*nwzy + jp*nwz + kp*nw + l ]/10000.0f*xWta) * yWta) * zWta) * wWtb +        //        f[ip][jp][kp][l ]*xWta*yWta*zWta*wWtb +
           ((((float)f[i *nwzy + j *nwz + k *nw + lp]/10000.0f*xWtb +                                //        f[i ][j ][k ][lp]*xWtb*yWtb*zWtb*wWta +
              (float)f[ip*nwzy + j *nwz + k *nw + lp]/10000.0f*xWta) * yWtb +                        //        f[ip][j ][k ][lp]*xWta*yWtb*zWtb*wWta +
             ((float)f[i *nwzy + jp*nwz + k *nw + lp]/10000.0f*xWtb +                                //        f[i ][jp][k ][lp]*xWtb*yWta*zWtb*wWta +
              (float)f[ip*nwzy + jp*nwz + k *nw + lp]/10000.0f*xWta) * yWta) * zWtb +                //        f[ip][jp][k ][lp]*xWta*yWta*zWtb*wWta +
            (((float)f[i *nwzy + j *nwz + kp*nw + lp]/10000.0f*xWtb +                                //        f[i ][j ][kp][lp]*xWtb*yWtb*zWta*wWta +
              (float)f[ip*nwzy + j *nwz + kp*nw + lp]/10000.0f*xWta) * yWtb +                        //        f[ip][j ][kp][lp]*xWta*yWtb*zWta*wWta +    
             ((float)f[i *nwzy + jp*nwz + kp*nw + lp]/10000.0f*xWtb +                                //        f[i ][jp][kp][lp]*xWtb*yWta*zWta*wWta +
              (float)f[ip*nwzy + jp*nwz + kp*nw + lp]/10000.0f*xWta) * yWta) * zWta) * wWta;         //        f[ip][jp][kp][lp]*xWta*yWta*zWta*wWta;
}

// Function to implement lower_bound
int lower_bound(float arr[], int N, float X)
{
    int mid;
 
    // Initialise starting index and
    // ending index
    int low = 0;
    int high = N;
 
    // Till low is less than high
    while (low < high) {
        mid = low + (high - low) / 2;
 
        // If X is less than or equal
        // to arr[mid], then find in
        // left subarray
        if (X <= arr[mid]) {
            high = mid;
        }
 
        // If X is greater arr[mid]
        // then find in right subarray
        else {
            low = mid + 1;
        }
    }
   
    // if X is greater than arr[n-1]
    if(low < N && arr[low] < X) {
       low++;
    }
       
    // Return the lower_bound index
    return low;
}
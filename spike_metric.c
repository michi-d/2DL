/*  Spike metric MEX-file (24.09.2012)
 *         written by Michael Drews
 *  
 *  Version log:
 *  24.09.2012: mex: birth
 *
 *  INPUT       : A, B, approx_steps, qt, qch, qa, c_add, c_del, one
 * 
 *
 *  A,B         : MATLAB arrays (max. 2-dimensional) which contain the amplitudes 
 *                of a spike signal at a certain time and in a certain channel.
 *                The time is given by the horizontal position in the array (second dimension: j)
 *                The channel is given by the vertical position in the array (first dimension: i)
 *  approx_steps: Level of approximation, number of time sections to build from the 
 *                arrays A and B. For approx_steps == 1 the result is exact.
 *  qt          : Cost of shifting an element by an index 1 in time direction. Cost 
 *                in relation to real time is then dependent from the sampling rate.
 *  qc          : Cost of shifting an element by an index 1 in channel direction. Cost
 *                in relation to frequency is then dependent from the frequency resolution.
 *  qa          : Cost of changing the amplitude of an element by a value of 1.
 *  c_add       : Cost of adding an element at an arbitrary position.
 *  c_del       : Cost of deleting an element at an arbitrary position (normally equals c_add).
 *  one         : Neutral element/value which indicates the absence of a spike signal (normally 0)
 *
 *
 *
 *  OUTPUT      : out
 *
 *  out         : An 1D-array which contains the distances of the corresponding time sections.
 *                Length is naturally given by approx_steps.
 *                For approx_steps == 1 the result is exact.
 *
 */



/*---------------------------------------------------------------------------*/
/*------------------------------------------- basic definitions -------------*/
/*---------------------------------------------------------------------------*/

#include "mex.h"
#include <math.h>


#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/*---------------------------------------------------------------------------*/
/*-------------------------------------- G L O B A L S ----------------------*/
/*---------------------------------------------------------------------------*/
     
double qt = 1.0;
double qc = 1.0;
double qa = 1.0;
double c_add = 1.0;
double c_del = 1.0;
double one   = 0.0;

/*---------------------------------------------------------------------------*/
/*------------------------------------------------ process network ----------*/
/*---------------------------------------------------------------------------*/


void print_array(double *array, int N, int M, int MAX_Z) {
    
 /* Prints out maximal 3-dimensional arrays with last dimension being interpreted
  * as third dimension.
  *
  * INPUT     : array, N, M, MAX_Z
  *
  * array     : array to be printed out
  * N         : integer containing the length of the first dimension of array
  * M         : integer containing the length of the second dimension of array
  * MAX_Z     : integer containing the length of the third dimension of array
  *
  * OUTPUT    : no output
  *
  */
  
  int i, j, z;
  double e;
  #define array(i,j,z) array[(i)+((j)+(z)*M)*N]
         
  for (z=0; z<MAX_Z; z++) {
      for (i=0; i<N; i++) {
          for (j=0; j<M; j++) {
             e = array(i,j,z);
             mexPrintf("%f ", e);
          }
          mexPrintf("\n");
      }
      mexPrintf("\n");
  }
  
}
     
    
double d_spike2D(double a1, double a2, double t1, double t2, double ch1, double ch2) {
 
 /*
  * INPUT     : a1, a2, t1, t2, ch1, ch2, mode
  *
  * a1        : amplitude of first element
  * a2        : amplitude of second element
  * t1        : time of first element
  * t2        : time of second element
  * ch1       : channel of first element
  * ch2       : channel of second element
  *
  * OUTPUT    : Returns the distance of two single elements to each other. 
  *
  */
    
    if ((a1 == one) && (a2 != one)) {
        return c_add;
    } else if ((a2 == one) && (a1 != one)) {
        return c_del;
    } else {
        return qa * fabs(a2 - a1) + qt * fabs(t2 - t1) + qc*fabs(ch2-ch1);
    }

}   
    
  
  void distance1D_slice(double *D_out, double *A, double *B, int N, int M, int U, int V, int fix1, int fix2, int mode) {
  
 /* Computes all possible distances of all possible subsequences of two 1-dimensional 
  * subsequences of A and B using the 1-dimensional algorithm.
  * The subsequences (or slices) are defined by fix1 and fix2. Depending on the mode
  * ((mode = 0) for vertical slices and (mode == 1) for horizontal slices) fix1 and fix2 define
  * either the vertical or horizontal slice positions in A and B.
  * All possible subsequences of these 2 slices are compared by saving the whole comparison matrix
  * D (given by the 1-dimensional algorithm) into D_out.
  *
  * INPUT     : D_out, A, B, N, M, U, V, fix1, fix2, mode
  *
  * D_out     : empty array of size NxMxUxV to store the comparison matrix (size NxU or MxV) in a part of it
  * A         : array containing the first spike signal
  * B         : array containing the second spike signal
  * N         : integer containing the length of the first dimension of A
  * M         : integer containing the length of the second dimension of A
  * U         : integer containing the length of the first dimension of B
  * V         : integer containing the length of the second dimension of B
  * fix1      : integer defining the slice position of the first subsequence
  *                  (fixed coordinate within A)
  * fix2      : integer defining the slice position of second the subsequence 
  *                  (fixed coordinate within B)
  * mode      : (mode == 0) for vertical slices and (mode == 1) for horizontal slices
  *
  *
  * OUTPUT    : no outputs / results are stored in D_out
  *
  */
      
  //mexPrintf("distance1D_slice...\n");   
  int i, j, k, l, N_tmp, M_tmp;
  double *A_tmp, *B_tmp, *D;
  
  if (mode == 0) {
    N_tmp = N;
    M_tmp = U;
  } else if (mode == 1) {
    N_tmp = M;
    M_tmp = V;
  }
    
  A_tmp = malloc( N_tmp*sizeof(double));
  B_tmp = malloc( M_tmp*sizeof(double));
  
  #define A(i,j) A[(i)+(j)*N]
  #define B(k,l) B[(k)+(l)*U]     
          
  if (mode == 0) {
      //mexPrintf("(%i, %i): \n", fix1, fix2);
      //mexPrintf("A: ", N_tmp);
      for (i=0; i<N_tmp; i++) {
         A_tmp[i] = A(i,fix1);
         //mexPrintf("%f, ", A_tmp[i]);
      }
      //mexPrintf("\n B: ");
      for (k=0; k<M_tmp; k++) {
         B_tmp[k] = B(k,fix2);
         //mexPrintf("%f, ", B_tmp[k]);
      }
      //mexPrintf("\n");
  } else if (mode == 1) {
      //mexPrintf("(%i, %i): \n", fix1, fix2);
      //mexPrintf("A: ", N_tmp);
      for (j=0; j<N_tmp; j++) {
         A_tmp[j] = A(fix1,j);
         //mexPrintf("%f, ", A_tmp[j]);
      }
      //mexPrintf("\n B: ");
      for (l=0; l<M_tmp; l++) {
         B_tmp[l] = B(fix2,l);
         //mexPrintf("%f, ", B_tmp[l]);
      }
      //mexPrintf("\n");
  }
  
  double *D1D;
  D1D = malloc(N_tmp*M_tmp*sizeof(double));
  
  #define D1D(i,j) D1D[(i)+(j)*N_tmp]
  #define D_out(i,j,k,l) D_out[(i)+((j)+((k)+(l)*U)*M)*N]
          
  double sum = 0;
  for (i=0; i<N_tmp; i++) {
    D1D(i,0) = sum + d_spike2D(A_tmp[i], one, 0.0, 0.0, 0.0, 0.0);
    sum = D1D(i,0);
    if (mode == 0) {
       D_out(i,fix1,0,fix2) = D1D(i,0);
    } else if (mode == 1) {
       D_out(fix1,i,fix2,0) = D1D(i,0);
    }
  }
  sum = 0;
  for (j=0; j<M_tmp; j++) {
    D1D(0,j) = sum + d_spike2D(one, B_tmp[j], 0.0, 0.0, 0.0, 0.0);
    sum = D1D(0,j);
    if (mode == 0) {
       D_out(0,fix1,j,fix2) = D1D(0,j);
    } else if (mode == 1) {
       D_out(fix1,0,fix2,j) = D1D(0,j);
    }
  }
  

  
  // run algorithm
 
  double t1, t2, ch1, ch2;

  for (i=1; i<N_tmp; i++) {
      for (j=1; j<M_tmp; j++) {
            if (mode == 0) {
               t1  = (double) fix1;
               t2  = (double) fix2;
               ch1 = (double) i;
               ch2 = (double) j;
            } else if (mode == 1) {
               t1  = (double) i;
               t2  = (double) j;
               ch1 = (double) fix1;
               ch2 = (double) fix2;
            }
  
          
            D1D(i,j) = min(D1D(i,j-1)    +  d_spike2D(one,      B_tmp[j], 0.0, 0.0, 0.0, 0.0),
                       min(D1D(i-1, j-1) +  d_spike2D(A_tmp[i], B_tmp[j], t1,  t2,  ch1, ch2),
                           D1D(i-1, j)   +  d_spike2D(A_tmp[i], one     , 0.0, 0.0, 0.0, 0.0) 
                          )
                          );
            
            D_out((int) ch1,(int) t1,(int) ch2,(int) t2) = D1D(i,j);
      }
  }
  
  //print_array(D1D, N_tmp, M_tmp, 1);  
  
  free(A_tmp);
  free(B_tmp);
  free(D1D);
    
}


void copy_pad(double *quelle, double *ziel, int Nq, int Mq, int N, int M, int end_j, int j_max, int s) {

 /* Copies an sub-array of quelle into an array ziel and paddes it with zeros on the left edges of 
  * every dimension.
  * Takes a time section out of quelle and prepares it for the 2-dimensional algorithm.
  *
  * INPUT     : quelle, ziel, Nq, Mq, N, M, end_j, j_max, s
  *
  * quelle    : array of size NqxMq to be copied
  * ziel      : empty array of size (Nq+1)x(Mq+1) 
  * Nq        : integer containing the length of the first dimension of quelle
  * Mq        : integer containing the length of the first dimension of quelle
  * N         : integer containing the length of the second dimension of ziel
  * M         : integer containing the length of the first dimension of ziel
  * end_j     : last index of current time section
  * j_max     : index width of one time sections
  * s         : current number of time section
  *
  *
  * OUTPUT    : no outputs / results are stored in ziel
  *
  */
    
  int i,j,z;
    
  //mexPrintf("copy_pad...\n");
  #define quelle(i,j,z) quelle[(i)+((j)+(z)*M)*N]
  #define ziel(i,j,z)   ziel[(i)+((j)+(z)*(Mq))*(Nq)]
    
  int MAX_Z = 1;
            
  for (z=0; z<MAX_Z; z++) {
     for (i=0; i<Nq; i++) {
         for (j=0; j<Mq; j++) {
             if ((i == 0) || (j == 0)) {
                 ziel(i,j,z) = one;
             } else {
                 ziel(i,j,z) = quelle(i - 1, j_max * s + (j-1), z);
             }
                 
        }
     }
  }
    
}


double zero_area(double *area, int N, int M, int start_i, int end_i, int start_j, int end_j, int pos) {
 
 /* Computes the cost of creating an sub-area of area (which means calculating 
  * the distance of every single element from the neutral element).
  *
  * INPUT     : area, N, M, start_i, end_i, start_j, end_j, pos
  *
  * area      : array of size NxM
  * N         : integer containing the length of the second dimension of area
  * M         : integer containing the length of the first dimension of area
  * start_i   : index defining the lower limit in the first dimension of the sub-area to be created
  * start_j   : index defining the lower limit in the second dimension of the sub-area to be created
  * end_i     : index defining the upper limit in the first dimension of the sub-area to be created
  * end_j     : index defining the upper limit in the second dimension of the sub-area to be created
  * pos       : position of sub-area in the metric funtion (d(a,0) or d(0,b)) which is the difference
  *             of assigning c_add or c_del as cost (creation of an element in B can be seen as a deletion in A)
  *             (pos == 1) means d(a,0) and (pos == 2) meand d(0,b)
  *
  * OUTPUT    : returns cost of creating the sub-area area(start_i ... end_i, start_j ... end_j).
  *
  */    
     
  double sum = 0;
  int i, j;
    
  //mexPrintf("d_zero_area...\n");
  #define area(i,j) area[(i)+(j)*N]
    
  for (i=start_i; i<end_i; i++) {
      for (j=start_j; j<end_j; j++) {
          if (pos == 1) {
              if (area(i, j) != one) {
                  sum += c_del;
              }
          } else if (pos == 2) {
              if (area(i, j) != one) {
                  sum += c_add;
              }
          }
      }
  }
  
  return sum;
}



void get_1D_distances(double *D_out, double *A, double*B, int N, int M,int U, int V, int mode) {
    
 /* Computes all one-dimensional distances which well be needed by the 2-dimensional algorithm.
  * Compares every possible horizontal slice of A with every possible horizontal slice of B and 
  * does the same for vertical slices too.
  *
  * INPUT     : D_out, A, B, N, M, U, V, mode
  *
  * D_out     : array of size NxMxUxV
  * A         : array containing the first spike signal
  * B         : array containing the second spike signal
  * N         : integer containing the length of the first dimension of A
  * M         : integer containing the length of the second dimension of A
  * U         : integer containing the length of the first dimension of B
  * V         : integer containing the length of the second dimension of B
  * mode      : (mode == 0) for vertical slices and (mode == 1) for horizontal slices
  *
  *
  * OUTPUT    : no output/results get stored in D_out laters
  *
  */    
    
    
  //mexPrintf("get_1D_distances...\n");
    
  if (mode == 0) {
      int j, l;
      for (j=0;j<M;j++) {
          for (l=0;l<V;l++) {
              distance1D_slice(D_out,A,B,N,M,U,V,j,l,mode);
          }
      }
  } else if (mode == 1) {
      int i, k;
      for (i=0;i<N;i++) {
          for (k=0;k<U;k++) {
              distance1D_slice(D_out,A,B,N,M,U,V,i,k,mode);
          }
      }
  }
    
}



void initialize(double *D, double *A, double *B, int N, int M, int U, int V) {
    
 /* Initializes the comparison matrix D used by the 2-dimensional algorithm.
  *
  *
  * INPUT     : D, A, B, N, M, U, V
  *
  * D         : empty array of size NxMxUxV, comparison matrix to be stored in
  * A         : array containing the first spike signal
  * B         : array containing the second spike signal
  * N         : integer containing the length of the first dimension of A
  * M         : integer containing the length of the second dimension of A
  * U         : integer containing the length of the first dimension of B
  * V         : integer containing the length of the second dimension of B
  *
  *
  * OUTPUT    : no output/results get stored in D
  *
  */    
    
  //mexPrintf("initialize...\n");
  int i,j,k,l;
  #define D(i,j,k,l) D[(i)+((j)+((k)+(l)*U)*M)*N]
            
  for (j=0; j<M; j++) {
      for (k=0; k<U; k++) {
          for (l=0;l<V;l++) {
              D(0,j,k,l) = zero_area(B, U, V, 1, k+1, 0, l+1, 2);
          }
      }
  }
    
  for (i=0; i<N; i++) {
      for (k=0; k<U; k++) {
          for (l=0;l<V;l++) {
              D(i,0,k,l) = zero_area(B, U, V, 0, k+1, 1, l+1, 2);
          }
      }
  }
    
  for (i=0; i<N; i++) {
      for (j=0; j<M; j++) {
          for (l=0;l<V;l++) {
              D(i,j,0,l) = zero_area(A, N, M, 1, i+1, 0, j+1, 1);
          }
      }
  }
    
  for (i=0; i<N; i++) {
      for (j=0; j<M; j++) {
          for (k=0;k<U;k++) {
              D(i,j,k,0) = zero_area(A, N, M, 0, i+1, 1, j+1, 1);
          }
      }
  }
    
       
}

double array_min(double *elements, int N) {
    
 /* Computes and returns the minimum of an array of N elements
  *
  *
  * INPUT     : elements, N
  *
  * elements  : array of N elements
  * N         : number of elements
  *
  *
  * OUTPUT    : value of the smallest element in elements
  *
  */    
    
  int n;
  double out = min(elements[0], elements[1]);
  for (n=2; n<N; n++) {
      out = min(out, elements[n]);
  }
  return out;
  
}


void get_costs(double *costs, double *D2D_, double *A, double *B, double *D1_, double *D2_, int i, int j, int k, int l, int N, int M, int U, int V) {

 /* Computes the cost of all possibilites to get to the element D(i,j,k,l) in the comparison
  * matrix D in the 2-dimensional algorithm. Returns an array of 15 elements because a cube
  * of length 2 has 16 units in 4 dimensions.
  *
  *
  * INPUT     : costs, D2D_, A, B, D1_, D2_, i, j, k, l, N, M, U, V
  *
  * costs     : array of 15 elements to store the calculated costs inside
  * D2D       : NxMxUxV array containing the comparison matrix in the 2-dimensional algorithm
  * A         : array containing the first spike signal
  * B         : array containing the second spike signal
  * i         : first dimension index of the element of D2D to be computed
  * j         : second dimension index of the element of D2D to be computed
  * k         : third dimension index of the element of D2D to be computed
  * l         : fourth dimension index of the element of D2D to be computed
  * N         : integer containing the length of the first dimension of A
  * M         : integer containing the length of the second dimension of A
  * U         : integer containing the length of the first dimension of B
  * V         : integer containing the length of the second dimension of B
  *
  *
  * OUTPUT    : no output/results get stored in costs
  *
  */   
    
  #define A(i,j) A[(i)+(j)*N]
  #define B(k,l) B[(k)+(l)*U]
  #define D1_(i,j,k,l) D1_[(i)+((j)+((k)+(l)*U)*M)*N]
  #define D2_(i,j,k,l) D2_[(i)+((j)+((k)+(l)*U)*M)*N]
  #define D2D_(i,j,k,l)   D2D_[(i)+((j)+((k)+(l)*U)*M)*N]
                         
  costs[0] = D2D_(i, j, k, l-1)  +  zero_area(B, U, V, 1, k+1, l, l+1, 2);
  costs[1] = D2D_(i, j-1, k, l)  +  zero_area(A, N, M, 1, i+1, j, j+1, 1);
  costs[2] = D2D_(i, j, k-1, l)  +  zero_area(B, U, V, k, k+1, 1, l+1, 2);
  costs[3] = D2D_(i-1, j, k, l)  +  zero_area(A, N, M, i, i+1, 1, j+1, 1);
    
  costs[4] = D2D_(i, j, k-1, l-1)  +  zero_area(B, U, V, k, k+1, 1, l+1, 2)  +  zero_area(B, U, V, 1, k-1+1, l, l+1, 2);
  costs[5] = D2D_(i-1, j-1, k, l)  +  zero_area(A, N, M, i, i+1, 1, j+1, 1)  +  zero_area(A, N, M, 1, i-1+1, j, j+1, 1);
    
  costs[6] = D2D_(i, j-1, k, l-1)  +  D1_(i,j,k,l);
  costs[7] = D2D_(i-1, j, k-1, l)  +  D2_(i,j,k,l);
    
  costs[8] = D2D_(i-1, j, k, l-1)  +  zero_area(B, U, V, 1, k+1, l, l+1, 2)  +  zero_area(A, N, M, i, i+1, 1, j+1, 1);
  costs[9] = D2D_(i, j-1, k-1, l)  +  zero_area(B, U, V, k, k+1, 1, l+1, 2)  +  zero_area(A, N, M, 1, i+1, j, j+1, 1);        
    
  costs[10]  =  D2D_(i, j-1, k-1, l-1)  +  zero_area(B, U, V, k, k+1, 1, l+1-1, 2)  +  D1_(i,j,k,l);
  costs[11]  =  D2D_(i-1, j, k-1, l-1)  +  zero_area(B, U, V, 1, k+1-1, l, l+1, 2)  +  D2_(i,j,k,l);
  costs[12]  =  D2D_(i-1, j-1, k, l-1)  +  zero_area(A, N, M, i, i+1, 1, j+1-1, 1)  +  D1_(i,j,k,l);
  costs[13]  =  D2D_(i-1, j-1, k-1, l)  +  zero_area(A, N, M, 1, i+1-1, j, j+1, 1)  +  D2_(i,j,k,l);
    
  costs[14]  =  D2D_(i-1, j-1, k-1, l-1)  +  min(D1_(i,j,k,l) + D2_(i,j-1,k,l-1), D1_(i-1,j,k-1,l) + D2_(i,j,k,l));
    
}


void distance2D(double *out, double *Ap, double *Bp, int N, int M, int U, int V, int S)
{
 /* Computes the distance of two (maximal 2-dimensional) array Ap and Bp via the Sellers' 
  * algorithm. Approximation level can be set by choosing S. For (S == 1) the algorithm is
  * exact.
  *
  *
  * INPUT     : out, Ap, Bp, N, M, U, B, S
  *
  * out       : empty array of size 1xS to store the results in
  * Ap        : array containing the first spike signal
  * Bb        : array containing the second spike signal
  * N         : integer containing the length of the first dimension of Ap
  * M         : integer containing the length of the second dimension of Ap
  * U         : integer containing the length of the first dimension of Bp
  * V         : integer containing the length of the second dimension of Bp
  * S         : integer indicating the approximation level to use / the number
  *                of seperate time sections to create from Ap and Bp
  *
  *
  * OUTPUT    : no output/results get stored in out
  *
  */   
    
  //mexPrintf("distance2D...\n");
 
  int i,j,k,l,z,s;
  double e;
 
  //print_array(Ap, N, M, MAX_Z);
 
  double *A_tmp, *B_tmp;
  double sum = 0.0;
  int j_max, l_max, end_j, end_l, N_tmp, M_tmp, U_tmp, V_tmp;
 
  j_max = (int)(M / S);
  l_max = (int)(V / S);
    
  for (s=0; s<S; s++) {
     end_j   = j_max * (s + 1);
     end_l   = l_max * (s + 1);
     if (s == S - 1) {
         end_j = M;
         end_l = V;
     }
    
     N_tmp = N+1;
     M_tmp = end_j - j_max*s + 1;
     U_tmp = U+1;
     V_tmp = end_l - l_max*s + 1;
    
     A_tmp = malloc( N_tmp*M_tmp*sizeof(double));
     B_tmp = malloc( U_tmp*V_tmp*sizeof(double));
    
     copy_pad(Ap, A_tmp, N_tmp, M_tmp, N, M, end_j, j_max, s);
     copy_pad(Bp, B_tmp, U_tmp, V_tmp, U, V, end_l, l_max, s);
    
     //print_array(A_tmp, N_tmp, M_tmp, 1);
     //print_array(B_tmp, U_tmp, V_tmp, 1);
    
     // initialize
    
     double *D2D;
     D2D = malloc( N_tmp*M_tmp*U_tmp*V_tmp*sizeof(double) );
     #define D2D(i,j,k,l) D2D[(i)+((j)+((k)+(l)*U_tmp)*M_tmp)*N_tmp]
    
     initialize(D2D, A_tmp, B_tmp, N_tmp, M_tmp, U_tmp, V_tmp);
         
     // get 1-D distances
              
     double *D1, *D2;
     D1 = malloc( N_tmp*M_tmp*U_tmp*V_tmp*sizeof(double) );
     D2 = malloc( N_tmp*M_tmp*U_tmp*V_tmp*sizeof(double) );
     #define D1(i,j,k,l) D1[(i)+((j)+((k)+(l)*U_tmp)*M_tmp)*N_tmp]
     #define D2(i,j,k,l) D2[(i)+((j)+((k)+(l)*U_tmp)*M_tmp)*N_tmp]
    
     get_1D_distances(D1, A_tmp, B_tmp, N_tmp, M_tmp, U_tmp, V_tmp, 0);
     get_1D_distances(D2, A_tmp, B_tmp, N_tmp, M_tmp, U_tmp, V_tmp, 1);
            
     // run algorithm
     double costs[15];
    
     for (i=1; i<N_tmp; i++) {
         for (j=1; j<M_tmp; j++) {
             for (k=1; k<U_tmp; k++) {
                 for (l=1; l<V_tmp; l++) {
                     get_costs(costs, D2D, A_tmp, B_tmp, D1, D2, i, j, k, l, N_tmp, M_tmp, U_tmp, V_tmp);
                     // mexPrintf("(%i, %i, %i, %i) : [%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f] \n", i, j, k, l, costs[0], costs[1], costs[2], costs[3], costs[4], costs[5], costs[6], costs[7], costs[8], costs[9], costs[10], costs[11], costs[12], costs[13], costs[14]);
                     D2D(i,j,k,l) = array_min(costs, 15);
                 }
             }
         }
     }
    
    
     //print D2D
     //for (i=0; i<N_tmp; i++) {
     //    for (j=0; j<M_tmp; j++) {
     //        for (k=0; k<U_tmp; k++) {
     //            for (l=0; l<V_tmp; l++) {
     //                mexPrintf("%f  ", D2D(i,j,k,l));
     //            }
     //            mexPrintf("\n");
     //        }
     //        mexPrintf("\n");
     //    }
     //   mexPrintf("\n");
     //}
    
     out[s] =  D2D(N_tmp - 1, M_tmp - 1, U_tmp - 1, V_tmp - 1);
     sum += out[s];
    
     free(A_tmp);
     free(B_tmp);
     free(D1);
     free(D2);
     free(D2D);
    
  }  
 
  // return sum;
         
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- M E X ----------------------------------*/
/*---------------------------------------------------------------------------*/


#include "mex.h"
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *Ap, *Bp, *dp;
  int     N, M, U, V, S, MAX_Z;
  
  /* Check for proper number of arguments. */
  if(nrhs!=9) {
    mexErrMsgTxt("9 inputs required: spike_time(:) = function(A(CHANNELS, SAMPLES), B(CHANNELS, SAMPLES), approx_level, qt, qch, qa, c_add, c_del, one)");
  } else if(nlhs>1) {
    mexErrMsgTxt("1 output: spike_time(:) = steps(approx_level)");
  }
  
  int *dims; 
  dims = (int*) mxGetDimensions(prhs[0]);
  N   = dims[0];
  M   = dims[1];

  dims = (int*) mxGetDimensions(prhs[1]);
  U   = dims[0];
  V   = dims[1];
  
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ) {
    mexErrMsgTxt("Input must be a noncomplex double matrix.");
  }
  
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
    mexErrMsgTxt("Input must be a noncomplex double matrix.");
  }
  
  /* Get scalar input f_s and weights*/
  S			= mxGetScalar(prhs[2]);
  qt        = mxGetScalar(prhs[3]);
  qc        = mxGetScalar(prhs[4]);
  qa        = mxGetScalar(prhs[5]);
  c_add     = mxGetScalar(prhs[6]);
  c_del     = mxGetScalar(prhs[7]);
  one       = mxGetScalar(prhs[8]);
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(S, 1, mxREAL);
  
  /* Assign pointers to each input and output. */
  Ap	= mxGetPr(prhs[0]);
  Bp	= mxGetPr(prhs[1]);
  dp	= mxGetPr(plhs[0]);

  distance2D(dp, Ap, Bp, N, M, U, V, S);

}
/*---------------------------------------------------------------------------*/
/*---------------------------------- E N D ----------------------------------*/
/*---------------------------------------------------------------------------*/

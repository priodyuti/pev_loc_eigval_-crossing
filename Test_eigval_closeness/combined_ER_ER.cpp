#include<iostream>
#include<fstream>
#include<cstdlib>
#include</home/priodyuti/Dropbox/RS/codes/Multiplex_Optimization/Header/dfs_template.h>
//#include</home/priodyuti/Dropbox/RS/codes/Multiplex_Optimization/Header/Eigen_template.h>
#include</home/priodyuti/Dropbox/RS/codes/Multiplex_Optimization/Header/random_template.h>
#include</home/priodyuti/Dropbox/RS/codes/Multiplex_Optimization/Header/eigen_all.h>
#include<math.h>
#include<ctime>
#include<stdlib.h>
#include<iomanip>
#include<sstream>
#include<string>
#include<algorithm>
#include<iterator>
#include<chrono>

using namespace std;
int main()
{
  unsigned int N1,N2,N,k;
 
  cout<<"Number of nodes in first ER comp: ";  cin>>N1; 
  cout<<"Number of nodes in second ER comp: "; cin >>N2;
  cout<<"Avg degree of the Network :"; cin>>k;
  N = N1 + N2;
  //unsigned int deg1 = (N1-1)*2;
  //unsigned int deg2 = N*k-deg1;
  //double k2 = deg2/double(N2);
  //cout<<"Deg of Star :"<<deg1<<' '<<"deg of ER :"<<deg2<<endl;
  //cout<<"Avg. degree of larger part :"<<k2<<endl;
  double p = k/double(N2);
  cout<<"Connection probability :"<<p<<endl; 
    
  bool** A1 = new bool*[N1];
  Allocate_2Darray(A1,N1,N1);
 
  bool** A2 = new bool*[N2];
  Allocate_2Darray(A2,N2,N2);
 
  bool** A = new bool*[N];
  Allocate_2Darray(A, N, N);
 
  unsigned long int Edges1 = ER_Random(N1,A1,p);
  //Display_Network(N1,A1);
 
  unsigned long int Edges2 = ER_Random(N2,A2,p);
   //Display_Network(N2,A2);
  //int N = N1 + N2;

  for(int i=0; i<N1; i++)
    for(int j=0; j<N1; j++) 
      A[i][j] = A1[i][j];

  A[N1-1][N1] = 1;
  A[N1][N1-1] = 1;
 //A[N1][N1+1]=1;
 //A[N1+1][N1]=1;

  for(int i=N1, m=0; i<N; i++,m++)
    for(int j=N1,n=0; j<N; j++,n++)
       A[i][j] = A2[m][n];  
 
  // Display_Network(N,A);
  double** Ev1 = new double*[N1];
  Allocate_2Darray(Ev1,N1,N1);
  double *Eval1 = new double[N1];

  Eigen_All(A1,N1,Ev1,Eval1);
 
  double c1 = IPR_Evec(Ev1,N1,N1-1);
  double c2 = IPR_Evec(Ev1,N1,N1-2);
  cout<<"IPR of x1 "<<c1<<' '<<"IPR of x2 :"<<c2<<endl;
  cout<<"lambda1: "<<Eval1[N1-1]<<' '<<"lambda2: "<<Eval1[N1-2]<<' '<<"lambda_3: "<<Eval1[N1-3]<<endl; 
  
  double** Ev2 = new double*[N2];
  Allocate_2Darray(Ev2,N2,N2);
  double *Eval2 = new double[N2];

  Eigen_All(A2,N2,Ev2,Eval2);
 
  c1 = IPR_Evec(Ev2,N2,N2-1);
  c2 = IPR_Evec(Ev2,N2,N2-2);
  cout<<"IPR of x1 "<<c1<<' '<<"IPR of x2 :"<<c2<<endl;
  cout<<"lambda_1: "<<Eval2[N2-1]<<' '<<"lambda_2: "<<Eval2[N2-2]<<' '<<"lambda_3: "<<Eval2[N2-3]<<endl; 
 
 //::::::::::::: Calculate IPR of A2 Matrix :::::::::::::::::
  double** Ev = new double*[N];
  Allocate_2Darray(Ev,N,N);
  double *Eval = new double[N];

  Eigen_All(A,N,Ev,Eval);
  /*cout<<endl;
  for(unsigned int i=0;i<N;i++)
    cout<<Eval[i]<<' ';
  cout<<endl<<endl;
  Display_Network(N,Ev);*/
 
  c1 = IPR_Evec(Ev,N,N-1);
  c2 = IPR_Evec(Ev,N,N-2);
  cout<<"IPR of x1 "<<c1<<' '<<"IPR of x2 :"<<c2<<endl;
  cout<<"lambda1: "<<Eval[N-1]<<' '<<"lambda2: "<<Eval[N-2]<<' '<<"lambda3: "<<Eval[N-3]<<endl; 

  unsigned int *Deg = new unsigned int[N];
  Degree_sequence(N,A,Deg);
  unsigned int max_deg = Max_Degree(A,N);
  unsigned int max_deg_index = Max_Degree_Index(A,N);
  unsigned int second_max_deg = Second_Max_Degree(A,N,max_deg_index,Deg);
  cout<<"max deg: "<<max_deg<<' '<<"second_max deg: "<<second_max_deg<<endl;
  
 //:::::::::::::Sorting Section for Eigen Vector::::::::::::::::::::::::::::::
  Network_File(A,N, "ER_ER_combined", 1);
 
  Free_2Darray(A,N);
  Free_2Darray(A1,N1);
  Free_2Darray(A2,N2);
  Free_2Darray(Ev,N);
  Free_2Darray(Ev1,N1);
  Free_2Darray(Ev2,N2);
  delete [] Eval;
  delete [] Eval1;
  delete [] Eval2;
 
 //cout<<"Number of Edges :"<<Edges<<endl;
 //cout<<"<k> :"<<(2*Edges)/double(N)<<endl;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 return 0;
}

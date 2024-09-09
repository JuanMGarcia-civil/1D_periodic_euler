#include <iostream> // std::cout
#include <fstream> // ofile
#include <string>
#include <math.h>
#include <iomanip> // set precision
#include <mpi.h>

#include "DataStructs.h"
#include "rk4.h"
#include "FluxFunctions.h"
#include "RHSoperator.h"

#ifdef _DOUBLE_
#define FLOATTYPE double
#else
#define FLOATTYPE float
#endif

// declare supporting functions
void write2File(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &U, std::string name);
FLOATTYPE calcL2norm(DataStruct<FLOATTYPE> &u, DataStruct<FLOATTYPE> &uinit);


int main(int narg, char **argv)
{
  
   int provided;
    MPI_Init_thread(&narg, &argv, MPI_THREAD_FUNNELED, &provided); // Inicialización de MPI
  
  int numPoints =  80;
  FLOATTYPE k = 2.; // wave number

  if(narg != 3) //nos dice que al ejecutar el programa hay que dar numero de puntos y numero ondas (en total tres argumentos)
  {
    std::cout<< "Wrong number of arguments. You should include:" << std::endl;
    std::cout<< "    Num points" << std::endl;
    std::cout<< "    Wave number" << std::endl;
    return 1;
  }else
  {
    numPoints = std::stoi(argv[1]);  // stoi string to integer
    k         = std::stod(argv[2]);  // string to double
  }

  // solution data
  DataStruct<FLOATTYPE> u(numPoints), f(numPoints), xj(numPoints);

  // flux function
  LinearFlux<FLOATTYPE> lf;

  // time solver
  RungeKutta4<FLOATTYPE> rk(u);

  // Initial Condition
  FLOATTYPE *datax = xj.getData();
  FLOATTYPE *dataU = u.getData();
  for(int j = 0; j < numPoints; j++)
  {
    // xj
    datax[j] = FLOATTYPE(j)/FLOATTYPE(numPoints-1);

    // init Uj
    dataU[j] = sin(k*2. * M_PI * datax[j]);
  }

  DataStruct<FLOATTYPE> Uinit;
  Uinit = u;

  // Operator
  Central1D<FLOATTYPE> rhs(u,xj,lf);

  FLOATTYPE CFL = 2.4;
  FLOATTYPE dt = CFL*datax[1];

  // Output Initial Condition
  write2File(xj, u, "initialCondition.csv");

  FLOATTYPE t_final = 1.; 
  FLOATTYPE time = 0.;
  DataStruct<FLOATTYPE> Ui(u.getSize()); // temp. data

  // init timer
  double compTime = MPI_Wtime();

  // main loop
  while(time < 1.) //quitada t_final, siempre es 1., optimiza ya que el compilador no evalua variable
  {
    if(time+dt >= 1.) dt = 1. - time;

    // take RK step
    rk.initRK();
  
    for(int s = 0; s < 4; s++) //quitada variable a leer del bucle
    {
      rk.stepUi(dt);
      Ui = *rk.currentU();
      rhs.eval(Ui);
      rk.setFi(rhs.ref2RHS());
    }
    rk.finalizeRK(dt);
    time += dt;
  }

  // finishe timer
  compTime = MPI_Wtime() - compTime;

  write2File(xj, u, "final.csv");

  // L2 norm
  FLOATTYPE err = calcL2norm(Uinit, u);
  std::cout << std::setprecision(4) << "Comp. time: " << compTime;
  std::cout << " sec. Error: " << err/k;
  std::cout << " kdx: " << k*datax[1]*2.*M_PI;
  std::cout << std::endl;

MPI_Finalize(); // Finalización de MPI

  return 0;
}


// ==================================================================
// AUXILIARY FUNCTIONS
// ==================================================================
void write2File(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &U, std::string name)
{
  std::ofstream file;
  file.open(name,std::ios_base::trunc);
  if(!file.is_open()) 
  {
    std::cout << "Couldn't open file for Initial Condition" << std::endl;
    exit(1);
  }
  
  for(int j = 0; j < U.getSize(); j++)
  {
    file << X.getData()[j] << " ," << U.getData()[j] << std::endl;
  }

  file.close();
}

FLOATTYPE calcL2norm(DataStruct<FLOATTYPE> &u, DataStruct<FLOATTYPE> &uinit)
{
  FLOATTYPE err = 0.;
  const FLOATTYPE *dataU = u.getData();
  const FLOATTYPE *dataInit = uinit.getData();

  for(int n = 0; n < u.getSize(); n++)
  {
    err += (dataU[n] - dataInit[n])*(dataU[n] - dataInit[n]);
  }

  return sqrt( err );
}

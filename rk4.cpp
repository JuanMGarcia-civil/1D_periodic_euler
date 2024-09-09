#include "rk4.h"
#include <cassert>

template<class T>
RungeKutta4<T>::RungeKutta4(DataStruct<T> &_Un):
Un(_Un)
{
  nSteps = 4;
  currentStep = 0;

  coeffsA = new T[4];
  coeffsB = new T[4];
  coeffsA[0] = 0.;
  coeffsA[1] = 0.5;
  coeffsA[2] = 0.5;
  coeffsA[3] = 1.;
  coeffsB[0] = 1.;
  coeffsB[1] = 2.;
  coeffsB[2] = 2.;
  coeffsB[3] = 1.;

  fi = new DataStruct<T>[2];  //dividimos fi para que tenga dimensión 2, variable menos pesada

  fi[0].setSize(Un.getSize());
  fi[1].setSize(Un.getSize());


  Ui.setSize(Un.getSize());
  Ufinal.setSize(Un.getSize()); // Para almacenar f1 + 2*f2 + 2*f3 + f4
};

template<class T>
RungeKutta4<T>::~RungeKutta4()
{
  delete[] fi;
  delete[] coeffsA;
  delete[] coeffsB;
};

template<class T>
int RungeKutta4<T>::getNumSteps()
{
  return nSteps;
};

template<class T>
void RungeKutta4<T>::initRK()
{
  currentStep = 0;

  // Inicializar Ufinal a cero
    T* dataUfinal = Ufinal.getData();
    for (int n = 0; n < Ufinal.getSize(); ++n) {
        dataUfinal[n] = 0.;
    }
};

template<class T>
void RungeKutta4<T>::stepUi(T dt)
{
  assert(currentStep < nSteps);

  if(currentStep == 0)
  {
    T *dataUi = Ui.getData();
    const T *dataU  = Un.getData();

    for(int n = 0; n < Ui.getSize(); n++)
    {
      dataUi[n] = dataU[n];
    }
  }
  else
  {
    T *datafi = fi[(currentStep-1)%2].getData(); //corregido al tamaño del nuevo array
    T *dataUi = Ui.getData();
    const T *dataU  = Un.getData();

    for(int n = 0; n < Ui.getSize(); n++)
    {
      dataUi[n] = dataU[n] + coeffsA[currentStep]*dt* datafi[n];
    }
  }
};

template<class T>
void RungeKutta4<T>::finalizeRK(const T dt)
{
  T *dataUn = Un.getData();
  const T* dataUfinal = Ufinal.getData(); //corregido

  const T oneDiv6 = 1. / 6.;
    for (int n = 0; n < Un.getSize(); ++n)
    {
        dataUn[n] += dt * oneDiv6 * dataUfinal[n];
    }
    
};

template<class T>
void RungeKutta4<T>::setFi(DataStruct<T> &_F)
{
    T* dataFi = fi[currentStep % 2].getData(); // Usando solo dos DataStruct
    const T* dataF = _F.getData();
    T* dataUfinal = Ufinal.getData();

    for (int n = 0; n < Ui.getSize(); ++n)
    {
        dataFi[n] = dataF[n];
        dataUfinal[n] += coeffsB[currentStep] * dataF[n];
    }

    currentStep++;
};

template<class T>
DataStruct<T> * RungeKutta4<T>::currentU()
{
  return &Ui;
};


template class RungeKutta4<float>;
template class RungeKutta4<double>;

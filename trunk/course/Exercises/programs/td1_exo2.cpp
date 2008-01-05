#include <iostream>
#include <fstream>
#include <blitz/array.h>

const int Dim = 2;

// Time steps
double delta_t = 0.1;

double spatialStep;

BZ_DECLARE_STENCIL2(heat_stencil,Tnp1,Tn)
Tnp1 += Tn/delta_t - Laplacian2D(Tn)/(spatialStep*spatialStep);
BZ_END_STENCIL





using namespace blitz;

typedef Array<double,Dim>   scalarField;

class HeatBC {
public:
  HeatBC()
  {}
  void apply(scalarField& T) const
  {
    // retrieve the dimension of the box
    int xN = T.ubound(firstDim);
    int yN = T.ubound(secondDim);
    Range all = Range::all();
    //apply Dirichlet BC on domain boundary
    // lower x
    T(0,all) = 0;
    // upper x
    T(xN,all) = 0;

    // lower y
    T(all,0) = 0;
    // upper y
    T(all,yN) = 0;
  }
};


template<typename T_stencil, typename T_numtype, int N_rank, typename T_BCs>
int conjugateGradientSolver(T_stencil stencil,
			    Array<T_numtype,N_rank>& x,
			    Array<T_numtype,N_rank>& rhs, double haltrho,
			    const T_BCs& boundaryConditions)
{
  // Calculate initial residual
  Array<T_numtype,N_rank> r = rhs.copy();
  r *= -1.0;

  boundaryConditions.apply(x);

  applyStencil(stencil, r, x);

  r *= -1.0;

  // Allocate the descent direction arrays
  Array<T_numtype,N_rank> p, q;
  allocateArrays(x.shape(), p, q);

  int iteration = 0;
  int converged = 0;
  T_numtype rho = 0.;
  T_numtype oldrho = 0.;

  const int maxIterations = 1000;

  while (iteration < maxIterations)
    {
      rho = sum(r * r);

      if ((iteration % 20) == 0)
	std::cout << "CG: Iter " << iteration << "\t rho = " << rho << std::endl;

      // Check halting condition
      if (rho < haltrho)
        {
	  converged = 1;
	  break;
        }

      if (iteration == 0)
        {
	  p = r;
        }
      else {
	T_numtype beta = rho / oldrho;
	p = beta * p + r;
      }
      q = 0.;
      //boundaryConditions.apply(p);
      applyStencil(stencil, q, p);
      T_numtype pq = sum(p*q);

      T_numtype alpha = rho / pq;
      x += alpha * p;
      r -= alpha * q;
      oldrho = rho;
      ++iteration;
    }

  if (!converged)
    std::cout << "Warning: CG solver did not converge" << std::endl;

  return iteration;
}



void snapshot(scalarField& T)
{
    static int snapshotNum = 0;

    ++snapshotNum;
    char filename[128];
    sprintf(filename, "temperature%03d.m", snapshotNum);

    ofstream ofs(filename);
    int N = T.length(firstDim);

    ofs << "T" << snapshotNum << " = [ ";
    for (int i=0; i < N; ++i)
    {
        for (int j=0; j < N; ++j)
        {
            float value = T(i,j);
            ofs << value << " ";
        }
        if (i < N-1)
            ofs << ";" << endl;
    }
    ofs << "];" << endl;
}

int main( int argc, char** argv )
{
  // Elapsed seconds
  double time_now = 0;


  // Arrays are NxN
  int N = 10;
  if ( argc == 2 )
    N = std::atoi( argv[1] );

  // A 1m x 1m x 1m domain
  spatialStep = 1.0 / (N - 1);



  // Allocate arrays: scalar fields
  scalarField T(N,N), nextT(N,N), F(N,N), rhs(N,N);

  T=0;
  F=0;
  F(Range(1,N-2),Range(1,N-2)) = 1;
  const int nIters = 1000;

  int i = 0;
  do
    {
      nextT = T;
      std::cout << "Iteration " << i << " Time = " << time_now << " s" 
		<< " ||T||_inf = " << max(abs(T)) 
		<< std::endl;

      rhs = F + T/delta_t;
      conjugateGradientSolver(heat_stencil(), T, rhs, 1e-8, HeatBC());

      double oldtime_now = time_now;
      time_now += delta_t;

      snapshot( T );
    }
  while ( sum((nextT-T)*(nextT-T)) > 1e-5);
}

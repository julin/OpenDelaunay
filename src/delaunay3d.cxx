#if defined(_OPENMP)
#include <omp.h>
#endif

#include <delaunay3d.h>

#include <stdio.h>
#include <stdlib.h>
#include <stack>
#include <algorithm>
#include <cstdint>

typedef double REAL;
extern void exactinit(int, int, int, REAL, REAL, REAL);
extern REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);
extern REAL insphere(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe);
extern REAL orient4d(REAL *pa, REAL *pb, REAL *pc, REAL *pd, REAL *pe,
  REAL ah, REAL bh, REAL ch, REAL dh, REAL eh);

DELAU_NS_BEGIN


#include "_private/delaunay_ds.h"

void BoundaryBox(const std::vector<double>& points, double bbox[2][3], double scale)
{
  bbox[0][0] = bbox[0][1] = bbox[0][2] = DBL_MAX;
  bbox[1][0] = bbox[1][1] = bbox[1][2] = -DBL_MAX;

  size_t nPoint = points.size() / 3;
  for (size_t i = 0; i < nPoint; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      bbox[0][j] = std::min(bbox[0][j], points[i * 3 + j]);
      bbox[1][j] = std::max(bbox[0][j], points[i * 3 + j]);
    }
  }

  if (fabs(scale - 1.0) > std::numeric_limits<double>::epsilon())
  {
    double factor = scale - 1.0;
    for (int j = 0; j < 3; j++)
    {
      double len = bbox[1][j] - bbox[0][j];
      len *= factor;
      bbox[0][j] -= len;
      bbox[1][j] += len;
    }
  }
}

void DelaunayTriangulation(const std::vector<double>& points, std::vector<size_t> &)
{
  enum
  {
    eInserted=0
  };
  std::vector<size_t> indices;
  SortHilbert(points, indices);
  double bbox[2][3];
  BoundaryBox(points, bbox, 1.01);
  exactinit(0, 1, 0, bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1], bbox[1][2] - bbox[0][2]);

  int nbBlocks = 1;
#if defined(_OPENMP)
  nbBlocks=omp_get_max_threads();
#endif

  std::vector<std::vector<size_t> > assignTo(nbBlocks);
  std::vector<size_t> assignTo0;
  std::vector<int> pointsColor(indices.size(),0);

  auto init = [&]() {
    size_t sizePerBlock = (nbBlocks * (indices.size() / nbBlocks)) / nbBlocks;
    size_t start = 0;
    for (int i = 0; i < nbBlocks; i++)
    {
      size_t end = start + sizePerBlock;
      assignTo[i].reserve(sizePerBlock);
      for (size_t j = start; j < end; j++)
      {
        assignTo[i].push_back(indices[j]);
        pointsColor[indices[j]] = i + 1;
      }
      start = end;
    }

    for (size_t j = start; j < indices.size(); j++)
    {
      assignTo[nbBlocks - 1].push_back(indices[j]);
      pointsColor[indices[j]] = nbBlocks;
    }

    srand((unsigned)points.size());
    size_t npoint0 = sizePerBlock /10;
    assignTo0.reserve(npoint0*nbBlocks);

    for (int i = 0; i < nbBlocks; i++)
    {
      for (size_t j = 0; j < npoint0; j++)
      {
         auto idx= rand() % assignTo[i].size();
         if (pointsColor[assignTo[i][idx]] != eInserted)
         {
           assignTo0.push_back(assignTo[i][idx]);
           pointsColor[assignTo[i][idx]] = eInserted;
         }
      }
    }
  };
  
  auto insertPoint = [&](const std::vector<size_t>& node2Insert) {

  };

  init();
  insertPoint(assignTo0);
}


DELAU_NS_END

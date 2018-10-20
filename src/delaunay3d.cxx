// Gmsh - Copyright (C) 1997-2018 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@onelab.info>.

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <delaunay3d.h>

#include <stdio.h>
#include <stdlib.h>
#include <stack>

DELAU_NS_BEGIN
void DelaunayTriangulation(const size_t numThreads, size_t arraysize, PointFunctor vertices, std::vector<size_t> &)
{
  std::vector<size_t> indices;
  SortHilbert(arraysize, vertices, indices);
  size_t nbBlocks = numThreads;

  std::vector<std::vector<size_t> > assignTo(nbBlocks);

  for (size_t i = 1; i < indices.size(); i++) {
    size_t start = indices[i - 1];
    size_t end = indices[i];
    size_t sizePerBlock = (nbBlocks * ((end - start) / nbBlocks)) / nbBlocks;
    // printf("sizePerBlock[%3d] = %8d\n",i,sizePerBlock);
    int currentBlock = 0;
    int localCounter = 0;
    {
      for (int jPt = start; jPt < end; jPt++) {
        if (localCounter++ >= sizePerBlock && currentBlock != nbBlocks - 1) {
          localCounter = 0;
          currentBlock++;
        }
        assignTo[currentBlock].push_back(jPt);
      }
    }
  }

}


DELAU_NS_END

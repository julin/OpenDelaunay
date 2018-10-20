#include <delaunay3d.h>
#include "spatial_sort.hxx"

#include <memory>
#include <algorithm>
#include <cmath>
#include <functional>

DELAU_NS_BEGIN
typedef std::function<void(size_t, double*)> PointFunctor;
void SortHilbert(size_t arraysize, PointFunctor vertices, std::vector<size_t> &indices);
void BoundaryBox(size_t arraysize, PointFunctor vertices, double bbox[2][3], double scale = 1.0);

void BoundaryBox(size_t nPoint, PointFunctor pointFun, double bbox[2][3], double scale)
{
  double pos[3];
  for (size_t i = 0; i < nPoint; i++) {
    pointFun(i, pos);
    for (int j = 0; j < 3; j++)
    {
      bbox[0][j] = std::min(bbox[0][j], pos[j]);
      bbox[1][j] = std::max(bbox[0][j], pos[j]);
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

struct HilbertSortB {
  // The code for generating table transgc from:
  // http://graphics.stanford.edu/~seander/bithacks.html.
  int transgc[8][3][8];
  int tsb1mod3[8];
  int maxDepth;
  int Limit;
  double bbox[2][3];
  PointFunctor pointFun;
  struct Vert
  {
    size_t pos;
  };

  void GetPosition(const Vert* vi,double* pos)const
  {
    pointFun(vi->pos,pos);
  }
  void ComputeGrayCode(int n);
  int Split(Vert **vertices, int arraysize, int GrayCode0, int GrayCode1,
    double BoundingBoxXmin, double BoundingBoxXmax,
    double BoundingBoxYmin, double BoundingBoxYmax,
    double BoundingBoxZmin, double BoundingBoxZmax);
  void Sort(Vert **vertices, int arraysize, int e, int d,
    double BoundingBoxXmin, double BoundingBoxXmax,
    double BoundingBoxYmin, double BoundingBoxYmax,
    double BoundingBoxZmin, double BoundingBoxZmax, int depth);
  HilbertSortB(PointFunctor vertices,int m = 0, int l = 2) : maxDepth(m), Limit(l)
  {
    pointFun = vertices;
    ComputeGrayCode(3);
  }
  void MultiscaleSortHilbert(Vert **vertices, int arraysize, int threshold,
    double ratio, int *depth,
    std::vector<size_t> &indices)
  {
    int middle = 0;
    if (arraysize >= threshold) {
      (*depth)++;
      middle = (int)(arraysize * ratio);
      MultiscaleSortHilbert(vertices, middle, threshold, ratio, depth, indices);
    }
    indices.push_back(middle);
    // printf("chunk starts at %d and size %d\n", middle, arraysize - middle);
    Sort(&(vertices[middle]), arraysize - middle, 0, 0, bbox[0][0],
      bbox[1][0], bbox[0][1], bbox[1][1], bbox[0][2],
      bbox[1][2], 0);
  }
  void Apply(std::vector<Vert *> &v, std::vector<size_t> &indices)
  {
    indices.clear();
    if (v.empty()) return;

    double pos[3];
    bbox[0][0] = bbox[0][1] = bbox[0][2] = DBL_MAX;
    bbox[1][0] = bbox[1][1] = bbox[1][2] = -DBL_MAX;
    for (size_t i = 0; i < v.size(); i++) {
      GetPosition(v[i],pos);
      for (int j = 0; j < 3; j++)
      {
        bbox[0][j] = std::min(bbox[0][j], pos[j]);
        bbox[1][j] = std::max(bbox[0][j], pos[j]);
      }
    }
    for (int j = 0; j < 3; j++)
    {
      double len = bbox[1][j] - bbox[0][j];
      len *= 0.01;
      bbox[0][j] -= len;
      bbox[1][j] += len;
    }

    Vert **pv = &v[0];
    int depth;
    indices.clear();
    MultiscaleSortHilbert(pv, (int)v.size(), 64, .125, &depth, indices);
    indices.push_back(v.size());
  }
};

void HilbertSortB::ComputeGrayCode(int n)
{
  int gc[8], N, mask, travel_bit;
  int e, d, f, k, g;
  int v, c;
  int i;

  N = (n == 2) ? 4 : 8;
  mask = (n == 2) ? 3 : 7;

  // Generate the Gray code sequence.
  for (i = 0; i < N; i++) {
    gc[i] = i ^ (i >> 1);
  }

  for (e = 0; e < N; e++) {
    for (d = 0; d < n; d++) {
      // Calculate the end point (f).
      f = e ^ (1 << d); // Toggle the d-th bit of 'e'.
                        // travel_bit = 2**p, the bit we want to travel.
      travel_bit = e ^ f;
      for (i = 0; i < N; i++) {
        // // Rotate gc[i] left by (p + 1) % n bits.
        k = gc[i] * (travel_bit * 2);
        g = ((k | (k / N)) & mask);
        // Calculate the permuted Gray code by xor with the start point (e).
        transgc[e][d][i] = (g ^ e);
      }
      //      assert(transgc[e][d][0] == e);
      //      assert(transgc[e][d][N - 1] == f);
    } // d
  } // e

    // Count the consecutive '1' bits (trailing) on the right.
  tsb1mod3[0] = 0;
  for (i = 1; i < N; i++) {
    v = ~i; // Count the 0s.
    v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
    for (c = 0; v; c++) {
      v >>= 1;
    }
    tsb1mod3[i] = c % n;
  }
}

int HilbertSortB::Split(Vert **vertices, int arraysize, int GrayCode0,
  int GrayCode1, double BoundingBoxXmin,
  double BoundingBoxXmax, double BoundingBoxYmin,
  double BoundingBoxYmax, double BoundingBoxZmin,
  double BoundingBoxZmax)
{
  Vert *swapvert;
  int axis, d;
  double split;

  // Find the current splitting axis. 'axis' is a value 0, or 1, or 2, which
  // correspoding to x-, or y- or z-axis.
  axis = (GrayCode0 ^ GrayCode1) >> 1;

  // Calulate the split position along the axis.
  if (axis == 0) {
    split = 0.5 * (BoundingBoxXmin + BoundingBoxXmax);
  }
  else if (axis == 1) {
    split = 0.5 * (BoundingBoxYmin + BoundingBoxYmax);
  }
  else { // == 2
    split = 0.5 * (BoundingBoxZmin + BoundingBoxZmax);
  }

  // Find the direction (+1 or -1) of the axis. If 'd' is +1, the direction of
  // the axis is to the positive of the axis, otherwise, it is -1.
  d = ((GrayCode0 & (1 << axis)) == 0) ? 1 : -1;

  // Partition the vertices into left- and right-arrays such that left points
  // have Hilbert indices lower than the right points.
  int i = 0;
  int j = arraysize - 1;
  double point[3];
  // Partition the vertices into left- and right-arrays.
  if (d > 0) {
    do {
      for (; i < arraysize; i++) {
        GetPosition(vertices[i], point);
        if (point[axis] >= split) {
          break;
        }
      }
      for (; j >= 0; j--) {
        GetPosition(vertices[j], point);
        if (point[axis] < split) 
          break;
      }
      // Is the partition finished?
      if (i == (j + 1)) 
        break;
      // Swap i-th and j-th vertices.
      swapvert = vertices[i];
      vertices[i] = vertices[j];
      vertices[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  }
  else {
    do {
      for (; i < arraysize; i++) {
        GetPosition(vertices[i], point);
        if (point[axis] <= split) 
          break;
      }
      for (; j >= 0; j--) {
        GetPosition(vertices[j], point);
        if (point[axis] > split) 
          break;
      }
      // Is the partition finished?
      if (i == (j + 1))
        break;
      // Swap i-th and j-th vertices.
      swapvert = vertices[i];
      vertices[i] = vertices[j];
      vertices[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  }

  return i;
}

// The sorting code is inspired by Tetgen 1.5
void HilbertSortB::Sort(Vert **vertices, int arraysize, int e, int d,
  double BoundingBoxXmin, double BoundingBoxXmax,
  double BoundingBoxYmin, double BoundingBoxYmax,
  double BoundingBoxZmin, double BoundingBoxZmax,
  int depth)
{
  double x1, x2, y1, y2, z1, z2;
  int p[9], w, e_w, d_w, k, ei, di;
  int n = 3, mask = 7;

  p[0] = 0;
  p[8] = arraysize;

  p[4] = Split(vertices, p[8], transgc[e][d][3], transgc[e][d][4],
    BoundingBoxXmin, BoundingBoxXmax, BoundingBoxYmin,
    BoundingBoxYmax, BoundingBoxZmin, BoundingBoxZmax);
  p[2] = Split(vertices, p[4], transgc[e][d][1], transgc[e][d][2],
    BoundingBoxXmin, BoundingBoxXmax, BoundingBoxYmin,
    BoundingBoxYmax, BoundingBoxZmin, BoundingBoxZmax);
  p[1] = Split(vertices, p[2], transgc[e][d][0], transgc[e][d][1],
    BoundingBoxXmin, BoundingBoxXmax, BoundingBoxYmin,
    BoundingBoxYmax, BoundingBoxZmin, BoundingBoxZmax);
  p[3] =
    Split(&(vertices[p[2]]), p[4] - p[2], transgc[e][d][2], transgc[e][d][3],
      BoundingBoxXmin, BoundingBoxXmax, BoundingBoxYmin, BoundingBoxYmax,
      BoundingBoxZmin, BoundingBoxZmax) +
    p[2];
  p[6] =
    Split(&(vertices[p[4]]), p[8] - p[4], transgc[e][d][5], transgc[e][d][6],
      BoundingBoxXmin, BoundingBoxXmax, BoundingBoxYmin, BoundingBoxYmax,
      BoundingBoxZmin, BoundingBoxZmax) +
    p[4];
  p[5] =
    Split(&(vertices[p[4]]), p[6] - p[4], transgc[e][d][4], transgc[e][d][5],
      BoundingBoxXmin, BoundingBoxXmax, BoundingBoxYmin, BoundingBoxYmax,
      BoundingBoxZmin, BoundingBoxZmax) +
    p[4];
  p[7] =
    Split(&(vertices[p[6]]), p[8] - p[6], transgc[e][d][6], transgc[e][d][7],
      BoundingBoxXmin, BoundingBoxXmax, BoundingBoxYmin, BoundingBoxYmax,
      BoundingBoxZmin, BoundingBoxZmax) +
    p[6];

  if (maxDepth > 0) {
    if ((depth + 1) == maxDepth) {
      // printf("max depth attained\n");
      return;
    }
  }

  // Recursively sort the points in sub-boxes.
  for (w = 0; w < 8; w++) {
    if ((p[w + 1] - p[w]) > Limit) {
      if (w == 0) {
        e_w = 0;
      }
      else {
        k = 2 * ((w - 1) / 2);
        e_w = k ^ (k >> 1);
      }
      k = e_w;
      e_w = ((k << (d + 1)) & mask) | ((k >> (n - d - 1)) & mask);
      ei = e ^ e_w;
      if (w == 0) {
        d_w = 0;
      }
      else {
        d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
      }
      di = (d + d_w + 1) % n;
      if (transgc[e][d][w] & 1) {
        x1 = 0.5 * (BoundingBoxXmin + BoundingBoxXmax);
        x2 = BoundingBoxXmax;
      }
      else {
        x1 = BoundingBoxXmin;
        x2 = 0.5 * (BoundingBoxXmin + BoundingBoxXmax);
      }
      if (transgc[e][d][w] & 2) { // y-axis
        y1 = 0.5 * (BoundingBoxYmin + BoundingBoxYmax);
        y2 = BoundingBoxYmax;
      }
      else {
        y1 = BoundingBoxYmin;
        y2 = 0.5 * (BoundingBoxYmin + BoundingBoxYmax);
      }
      if (transgc[e][d][w] & 4) { // z-axis
        z1 = 0.5 * (BoundingBoxZmin + BoundingBoxZmax);
        z2 = BoundingBoxZmax;
      }
      else {
        z1 = BoundingBoxZmin;
        z2 = 0.5 * (BoundingBoxZmin + BoundingBoxZmax);
      }
      Sort(&(vertices[p[w]]), p[w + 1] - p[w], ei, di, x1, x2, y1, y2, z1, z2,
        depth + 1);
    }
  }
}

void SortHilbert(const std::vector<double>& points, std::vector<size_t> &indices)
{
  class SpatialSortingTraits
  {
    const std::vector<double>& pointFun;
  public:
    SpatialSortingTraits(const std::vector<double>& ptFun)
      :pointFun(ptFun) {}

    typedef size_t Point_3;

    struct Less_x_3
    {
      const std::vector<double>& pointFun;

      Less_x_3(const std::vector<double>& ptFun)
        :pointFun(ptFun) {}

      bool operator()(size_t p, size_t q) const
      {
        return pointFun[p*3] < pointFun[q * 3];
      }
    };

    struct Less_y_3
    {
      const std::vector<double>& pointFun;
      Less_y_3(const std::vector<double>& ptFun)
        :pointFun(ptFun) {}

      bool operator()(size_t p, size_t q) const
      {
        return pointFun[p * 3 + 1] < pointFun[q * 3 + 1];
      }
    };

    struct Less_z_3
    {
      const std::vector<double>& pointFun;
      Less_z_3(const std::vector<double>& ptFun)
        :pointFun(ptFun) {}

      bool operator()(size_t p, size_t q) const
      {
        return pointFun[p * 3 + 2] < pointFun[q * 3 + 2];
      }
    };


    Less_x_3 less_x_3_object() const { return Less_x_3(pointFun); }

    Less_y_3 less_y_3_object() const { return Less_y_3(pointFun); }

    Less_z_3 less_z_3_object() const { return Less_z_3(pointFun); }
  };
  indices.resize(points.size()/3);
  for (size_t i = 0; i < indices.size(); i++)
  {
    indices[i] = i;
  }
  SpatialSortingTraits gt(points);
  spatial_sort(indices.begin(), indices.end(), gt);
}

void SortHilbert(size_t arraysize, PointFunctor vertices, std::vector<size_t> &indices)
{
  class SpatialSortingTraits
  {
    PointFunctor pointFun;
  public:
    SpatialSortingTraits(PointFunctor ptFun)
      :pointFun(ptFun){}

    typedef size_t Point_3;

    struct Less_x_3
    {
      PointFunctor pointFun;
     
      Less_x_3(PointFunctor ptFun)
        :pointFun(ptFun) {}

      bool operator()(size_t p, size_t q) const
      {
        double x[3], y[3];
        pointFun(p, x);
        pointFun(q, y);
        return x[0] < y[0];
      }
    };

    struct Less_y_3
    {
      PointFunctor pointFun;
      Less_y_3(PointFunctor ptFun)
        :pointFun(ptFun) {}

      bool operator()(size_t p, size_t q) const
      {
        double x[3], y[3];
        pointFun(p, x);
        pointFun(q, y);
        return x[1] < y[1];
      }
    };

    struct Less_z_3
    {
      PointFunctor pointFun;
      Less_z_3(PointFunctor ptFun)
        :pointFun(ptFun) {}

      bool operator()(size_t p, size_t q) const
      {
        double x[3], y[3];
        pointFun(p, x);
        pointFun(q, y);
        return x[2] < y[2];
      }
    };


    Less_x_3 less_x_3_object() const { return Less_x_3(pointFun); }

    Less_y_3 less_y_3_object() const { return Less_y_3(pointFun); }

    Less_z_3 less_z_3_object() const { return Less_z_3(pointFun); }
  };

  if (0 ) {
    indices.resize(arraysize);
    for (size_t i = 0; i < arraysize; i++)
    {
      indices[i] = i;
    }
    SpatialSortingTraits gt(vertices);
    spatial_sort(indices.begin(), indices.end(), gt);
  }
  else {
    std::vector<HilbertSortB::Vert> vReal(arraysize);
    std::vector<HilbertSortB::Vert *> v(arraysize);
    HilbertSortB h(vertices,1000);
    for (size_t i = 0; i < arraysize; i++)
    {
      vReal[i].pos = i;
      v[i] = &(vReal[i]);
    }
    h.Apply(v, indices);
  }
}

DELAU_NS_END


//&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined(_ENABLE_UNITEST)
#include "test/test_config.h"
TEST(SortHilbertTest, test_SortHilbert)
{
  MEXT_ENABLE_ASSERT;
  unsigned nPoint = 20000000;
  srand(nPoint);
  std::vector<size_t> indices;
  std::vector<double> points(nPoint * 3);
  for (size_t i=0;i<points.size();i++)
  {
    points[i] = rand() % 5000;
  }

  DELAU_NS::SortHilbert(points, indices);

  int xx = 0;
}
#endif
 
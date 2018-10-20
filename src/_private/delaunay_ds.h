#pragma once
#ifndef MAX_NUM_THREADS_
#define MAX_NUM_THREADS_ 8
#endif

template <class T> class PageAlloc {
public:
  std::vector<T *> _all;
  size_t _current;
  size_t _nbAlloc;
  size_t size() { return _current + (_all.size() - 1) * _nbAlloc; }

public:
  T * operator()(size_t i)
  {
    const size_t _array = i / _nbAlloc;
    const size_t _offset = i % _nbAlloc;
    // printf("%d %d %d\n",i,_array,_offset);
    return _all[_array] + _offset;
  }

  PageAlloc(unsigned int s) : _current(0), _nbAlloc(s)
  {
    _all.push_back(new T[_nbAlloc]);
  }

  ~PageAlloc()
  {
    for (unsigned int i = 0; i < _all.size(); i++) {
      delete[] _all[i];
    }
  }

  T *newStuff()
  {
    if (_current == _nbAlloc) {
      _all.push_back(new T[_nbAlloc]);
      _current = 0;
    }
    _current++;
    return _all[_all.size() - 1] + (_current - 1);
  }
};

typedef uint8_t CHECKTYPE;
struct Tetra;

struct Vert
{
private:
  const double* _x;
  double _lc;
  size_t _num;
  Tetra   *_t;
public:
  void setT(Tetra *t) { _t = t; }
  Tetra *getT() const { return _t; }
  size_t getNum() const { return _num; }
  void setNum(unsigned int n) { _num = n; }

  double operator[](int i)const { return _x[i]; }
  double lc() const { return _lc; }
  double &lc() { return _lc; }
  inline operator const double *() { return _x; }

  Vert(const double* x = nullptr, double lc = 0, int num = 0)
    :_x(x),_lc(lc), _num(num), _t(nullptr)
  {
    thread = 0;
  }
  uint8_t thread;
};

inline bool inSphereTest_s(Vert *va, Vert *vb, Vert *vc, Vert *vd, Vert *ve)
{
  double val = insphere(
    (double *)va, (double *)vb, (double *)vc, (double *)vd, (double *)ve);
  if (val == 0.0) {
    int count;
    // symbolic perturbation
    Vert *pt[5] = { va, vb, vc, vd, ve };
    int swaps = 0;
    int n = 5;
    do {
      count = 0;
      n = n - 1;
      for (int i = 0; i < n; i++) {
        if (pt[i] > pt[i + 1]) {
          Vert *swappt = pt[i];
          pt[i] = pt[i + 1];
          pt[i + 1] = swappt;
          count++;
        }
      }
      swaps += count;
    } while (count > 0);
    double oriA = orient3d((double *)pt[1], (double *)pt[2],
      (double *)pt[3], (double *)pt[4]);
    if (oriA != 0.0) {
      // Flip the sign if there are odd number of swaps.
      if ((swaps % 2) != 0) oriA = -oriA;
      val = oriA;
    }
    else {
      double oriB = -orient3d(
        (double *)pt[0], (double *)pt[2], (double *)pt[3], (double *)pt[4]);
      // Flip the sign if there are odd number of swaps.
      if ((swaps % 2) != 0) oriB = -oriB;
      val = oriB;
    }
  }
  return val > 0;
}

class Edge {
public:
  Vert * first, *second;
  Edge(Vert *v1, Vert *v2) : first(std::min(v1, v2)), second(std::max(v1, v2))
  {
  }
  bool operator==(const Edge &e) const
  {
    return e.first == first && e.second == second;
  }
  bool operator<(const Edge &e) const
  {
    if (first < e.first) return true;
    if (first > e.first) return false;
    if (second < e.second) return true;
    return false;
  }
};

struct EdgeContainer {
public:
  std::vector<std::vector<Edge> > _hash;
  std::size_t _size, _size_obj;

public:
  EdgeContainer(size_t N = 1000000)
  {
    _size = 0;
    _hash.resize(N);
    _size_obj = sizeof(Edge);
  }

  std::size_t H(const Edge &edge) const
  {
    const std::size_t h = ((std::size_t)edge.first);
    //    printf("%lu %lu %lu %lu\n",h,(h/2)%_hash.size(),h/64,h>>6);
    return (h / _size_obj) % _hash.size();
  }

  bool find(const Edge &e) const
  {
    std::vector<Edge> const &v = _hash[H(e)];
    return std::find(v.begin(), v.end(), e) != v.end();
  }

  bool empty() const { return _size == 0; }

  bool addNewEdge(const Edge &e)
  {
    std::vector<Edge> &v = _hash[H(e)];
    for (unsigned int i = 0; i < v.size(); i++)
      if (e == v[i]) {
        return false;
      }
    v.push_back(e);
    _size++;
    return true;
  }
};


template<typename T>
void cswap(T& a,T& b)
{
  if (a > b) {
    std::swap(a,b);
  }
}


struct Face {
  Vert *v[3];
  Vert *V[3];
  Face(Vert *v1, Vert *v2, Vert *v3)
  {
    V[0] = v[0] = v1;
    V[1] = v[1] = v2;
    V[2] = v[2] = v3;
    cswap(v[0], v[1]);
    cswap(v[1], v[2]);
    cswap(v[0], v[1]);
  }

  bool operator==(const Face &other) const
  {
    return v[0] == other.v[0] && v[1] == other.v[1] && v[2] == other.v[2];
  }

  bool operator<(const Face &other) const
  {
    if (v[0] < other.v[0]) return true;
    if (v[0] > other.v[0]) return false;
    if (v[1] < other.v[1]) return true;
    if (v[1] > other.v[1]) return false;
    if (v[2] < other.v[2]) return true;
    return false;
  }
};


struct Tetra 
{
  Tetra *T[4];
  Vert *V[4];
  CHECKTYPE _bitset[MAX_NUM_THREADS_];
  bool _modified;
  //  static int in_sphere_counter;
  Tetra() : _modified(true)
  {
    V[0] = V[1] = V[2] = V[3] = nullptr;
    T[0] = T[1] = T[2] = T[3] = nullptr;
    setAllDeleted();
  }
  //  inline bool isFace () const {return V[3]==NULL;}
  int setVerticesNoTest(Vert *v0, Vert *v1, Vert *v2, Vert *v3)
  {
    _modified = true;
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
    for (int i = 0; i < 4; i++) {
      if (V[i]) V[i]->setT(this);
    }
    //    for (int i=0;i<4;i++)_copy[i] = *V[i];
    return 1;
  }
  int setVertices(Vert *v0, Vert *v1, Vert *v2, Vert *v3)
  {
    _modified = true;
    double val = orient3d((double *)v0, (double *)v1,(double *)v2, (double *)v3);
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
    for (int i = 0; i < 4; i++) {
      if (V[i]) 
        V[i]->setT(this);
    }
    if (val > 0) {
      // for (int i=0;i<4;i++)_copy[i] = *V[i];
      return 1;
    }
    else if (val < 0) {
      // throw;
      V[0] = v1;
      V[1] = v0;
      V[2] = v2;
      V[3] = v3;
      // for (int i=0;i<4;i++)_copy[i] = *V[i];
      return -1;
    }
    else {
      // throw;
      return 0;
    }
  }

  Tetra(Vert *v0, Vert *v1, Vert *v2, Vert *v3)
  {
    setVertices(v0, v1, v2, v3);
    T[0] = T[1] = T[2] = T[3] = nullptr;
    setAllDeleted();
  }

  void setAllDeleted()
  {
    for (int i = 0; i < MAX_NUM_THREADS_; i++) 
      _bitset[i] = 0x0;
  }

  void unset(int thread, int iPnt) { _bitset[thread] &= ~(1 << iPnt); }

  void set(int thread, int iPnt) { _bitset[thread] |= (1 << iPnt); }

  CHECKTYPE isSet(int thread, int iPnt) const
  {
    return _bitset[thread] & (1 << iPnt);
  }

  Face getFace(int k) const
  {
    const int fac[4][3] = { { 0, 1, 2 },{ 1, 3, 2 },{ 2, 3, 0 },{ 1, 0, 3 } };
    return Face(V[fac[k][0]], V[fac[k][1]], V[fac[k][2]]);
  }

  Vert *getOppositeVertex(int k) const
  {
    const int o[4] = { 3, 0, 1, 2 };
    return V[o[k]];
  }

  Edge getEdge(int k) const
  {
    const int edg[6][2] = { { 0, 1 },{ 0, 2 },{ 0, 3 },{ 1, 2 },{ 1, 3 },{ 2, 3 } };
    return Edge(std::min(V[edg[k][0]], V[edg[k][1]]),
      std::max(V[edg[k][0]], V[edg[k][1]]));
  }

  bool inSphere(Vert *vd, int thread)
  {
    //    in_sphere_counter++;
    return inSphereTest_s(V[0], V[1], V[2], V[3], vd);
  }
};

struct Conn {
  Face f;
  int i;
  Tetra *t;
  Conn() : f(0, 0, 0), i(0), t(0) {}
  Conn(Face _f, int _i, Tetra *_t) : f(_f), i(_i), t(_t) {}
  bool operator==(const Conn &c) const { return f == c.f; }
  bool operator<(const Conn &c) const { return f < c.f; }
};

class TetraContainer {
  std::vector<PageAlloc<Tetra> *> _perThread;

public:
  size_t size(int thread) const
  {
    if ((int)_perThread.size() <= thread) return 0;
    return _perThread[thread]->size();
  }

  Tetra *operator()(int thread, int j) const { return (*_perThread[thread])(j); }

  TetraContainer(int nbThreads, int preallocSizePerThread)
  {
    _perThread.resize(nbThreads);
#if defined(_OPENMP)
#pragma omp parallel num_threads(nbThreads)
#endif
    {
#if defined(_OPENMP)
      int myThread = omp_get_thread_num();
#else
      int myThread = 0;
#endif
      _perThread[myThread] = new PageAlloc<Tetra>(preallocSizePerThread);
    }
  }
  Tetra *newTet(int thread) { return _perThread[thread]->newStuff(); }
  ~TetraContainer() { delete _perThread[0]; }
};

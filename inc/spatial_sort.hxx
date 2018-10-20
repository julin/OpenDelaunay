#pragma once
#include <algorithm>
#include <iterator>
#include <cstddef>
#include <cassert>

namespace internal {
    template <class K, int x, bool up> struct Hilbert_cmp_3;

    template <class K, int x>
    struct Hilbert_cmp_3<K,x,true>
      : public std::binary_function<typename K::Point_3,
      typename K::Point_3, bool>
    {
      typedef typename K::Point_3 Point;
      K k;
      Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
      bool operator() (const Point &p, const Point &q) const
      {
        return Hilbert_cmp_3<K,x,false> (k) (q, p);
      }
    };

    template <class K>
    struct Hilbert_cmp_3<K,0,false>
      : public std::binary_function<typename K::Point_3,
      typename K::Point_3, bool>
    {
      typedef typename K::Point_3 Point;
      K k;
      Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
      bool operator() (const Point &p, const Point &q) const
      {
        return k.less_x_3_object() (p, q);
      }
    };

    template <class K>
    struct Hilbert_cmp_3<K,1,false>
      : public std::binary_function<typename K::Point_3,
      typename K::Point_3, bool>
    {
      typedef typename K::Point_3 Point;
      K k;
      Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
      bool operator() (const Point &p, const Point &q) const
      {
        return k.less_y_3_object() (p, q);
      }
    };

    template <class K>
    struct Hilbert_cmp_3<K,2,false>
      : public std::binary_function<typename K::Point_3,
      typename K::Point_3, bool>
    {
      typedef typename K::Point_3 Point;
      K k;
      Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
      bool operator() (const Point &p, const Point &q) const
      {
        return k.less_z_3_object() (p, q);
      }
    };

    template <class RandomAccessIterator, class Cmp>
    RandomAccessIterator
    hilbert_split (RandomAccessIterator begin, RandomAccessIterator end,
                   Cmp cmp = Cmp ())
    {
        if (begin >= end) return begin;

        RandomAccessIterator middle = begin + (end - begin) / 2;
        std::nth_element (begin, middle, end, cmp);
        return middle;
    }
}

template <class K>
class Hilbert_sort_median_3
{
public:
  typedef K Kernel;
  typedef typename Kernel::Point_3 Point;

private:
  Kernel _k;
  std::ptrdiff_t _limit;

  template <int x, bool up> struct Cmp : public internal::Hilbert_cmp_3<Kernel,x,up>
  { Cmp (const Kernel &k) : internal::Hilbert_cmp_3<Kernel,x,up> (k) {} };

public:
  Hilbert_sort_median_3 (const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
    : _k(k), _limit (limit)
  {}

  template <int x, bool upx, bool upy, bool upz, class RandomAccessIterator>
  void sort (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    const int y = (x + 1) % 3, z = (x + 2) % 3;
    if (end - begin <= _limit) return;

    RandomAccessIterator m0 = begin, m8 = end;

    RandomAccessIterator m4 = internal::hilbert_split (m0, m8, Cmp< x,  upx> (_k));
    RandomAccessIterator m2 = internal::hilbert_split (m0, m4, Cmp< y,  upy> (_k));
    RandomAccessIterator m1 = internal::hilbert_split (m0, m2, Cmp< z,  upz> (_k));
    RandomAccessIterator m3 = internal::hilbert_split (m2, m4, Cmp< z, !upz> (_k));
    RandomAccessIterator m6 = internal::hilbert_split (m4, m8, Cmp< y, !upy> (_k));
    RandomAccessIterator m5 = internal::hilbert_split (m4, m6, Cmp< z,  upz> (_k));
    RandomAccessIterator m7 = internal::hilbert_split (m6, m8, Cmp< z, !upz> (_k));

    sort<z, upz, upx, upy> (m0, m1);
    sort<y, upy, upz, upx> (m1, m2);
    sort<y, upy, upz, upx> (m2, m3);
    sort<x, upx,!upy,!upz> (m3, m4);
    sort<x, upx,!upy,!upz> (m4, m5);
    sort<y,!upy, upz,!upx> (m5, m6);
    sort<y,!upy, upz,!upx> (m6, m7);
    sort<z,!upz,!upx, upy> (m7, m8);
  }

  template <class RandomAccessIterator>
  void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    sort <0, false, false, false> (begin, end);
  }
};

namespace internal {
  template <class K, int x, bool up> struct Fixed_hilbert_cmp_3;

  template <class K, int x>
  struct Fixed_hilbert_cmp_3<K,x,true>
    : public std::binary_function<typename K::Point_3,
    typename K::Point_3, bool>
  {
    typedef typename K::Point_3 Point;
    K k;
    double value;
    Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
    bool operator() (const Point &p) const
    {
      return ! Fixed_hilbert_cmp_3<K,x,false> (value,k) (p);
    }
  };

  template <class K>
  struct Fixed_hilbert_cmp_3<K,0,false>
    : public std::binary_function<typename K::Point_3,
    typename K::Point_3, bool>
  {
    typedef typename K::Point_3 Point;
    K k;
    double value;
    Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
    bool operator() (const Point &p) const
    {
      return to_double(k.compute_x_3_object()(p)) < value;
    }
  };

  template <class K>
  struct Fixed_hilbert_cmp_3<K,1,false>
    : public std::binary_function<typename K::Point_3,
    typename K::Point_3, bool>
  {
    typedef typename K::Point_3 Point;
    K k;
    double value;
    Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
    bool operator() (const Point &p) const
    {
      return to_double(k.compute_y_3_object()(p)) < value;
    }
  };

  template <class K>
  struct Fixed_hilbert_cmp_3<K,2,false>
    : public std::binary_function<typename K::Point_3,
    typename K::Point_3, bool>
  {
    typedef typename K::Point_3 Point;
    K k;
    double value;
    Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
    bool operator() (const Point &p) const
    {
      return to_double(k.compute_z_3_object()(p)) < value ;
    }
  };

  template <class RandomAccessIterator, class Cmp>
  RandomAccessIterator
    fixed_hilbert_split (RandomAccessIterator begin, RandomAccessIterator end,
    Cmp cmp = Cmp ())
  {
    if (begin >= end) return begin;

    return std::partition (begin, end, cmp);
  }
}

template <class K>
class Hilbert_sort_middle_3
{
public:
  typedef K Kernel;
  typedef typename Kernel::Point_3 Point;

private:
  Kernel _k;
  std::ptrdiff_t _limit;

  template <int x, bool up> struct Cmp : public internal::Fixed_hilbert_cmp_3<Kernel,x,up>
  { Cmp (double v,const Kernel &k) : internal::Fixed_hilbert_cmp_3<Kernel,x,up> (v,k) {} };

public:
  Hilbert_sort_middle_3 (const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
    : _k(k), _limit (limit)
  {}

  template <int x, bool upx, bool upy, bool upz, class RandomAccessIterator>
  void sort (RandomAccessIterator begin, RandomAccessIterator end,
    double xmin, double ymin, double zmin, 
    double xmax, double ymax, double zmax) const
  {
    const int y = (x + 1) % 3, z = (x + 2) % 3;
    if (end - begin <= _limit) return;

    double xmed= (xmin+xmax)/2;
    double ymed= (ymin+ymax)/2;
    double zmed= (zmin+zmax)/2;


    RandomAccessIterator m0 = begin, m8 = end;

    RandomAccessIterator m4 = 
      internal::fixed_hilbert_split (m0, m8, Cmp< x,  upx> (xmed,_k));
    RandomAccessIterator m2 = 
      internal::fixed_hilbert_split (m0, m4, Cmp< y,  upy> (ymed,_k));
    RandomAccessIterator m6 = 
      internal::fixed_hilbert_split (m4, m8, Cmp< y, !upy> (ymed,_k));
    RandomAccessIterator m1 = 
      internal::fixed_hilbert_split (m0, m2, Cmp< z,  upz> (zmed,_k));
    RandomAccessIterator m3 = 
      internal::fixed_hilbert_split (m2, m4, Cmp< z, !upz> (zmed,_k));
    RandomAccessIterator m5 = 
      internal::fixed_hilbert_split (m4, m6, Cmp< z,  upz> (zmed,_k));
    RandomAccessIterator m7 = 
      internal::fixed_hilbert_split (m6, m8, Cmp< z, !upz> (zmed,_k));


    sort<z, upz, upx, upy> (m0, m1, zmin, xmin, ymin, zmed, xmed, ymed);
    sort<y, upy, upz, upx> (m1, m2, ymin, zmed, xmin, ymed, zmax, xmed);
    sort<y, upy, upz, upx> (m2, m3, ymed, zmed, xmin, ymax, zmax, xmed);
    sort<x, upx,!upy,!upz> (m3, m4, xmin, ymax, zmed, xmed, ymed, zmin);
    sort<x, upx,!upy,!upz> (m4, m5, xmed, ymax, zmed, xmax, ymed, zmin);
    sort<y,!upy, upz,!upx> (m5, m6, ymax, zmed, xmax, ymed, zmax, xmed);
    sort<y,!upy, upz,!upx> (m6, m7, ymed, zmed, xmax, ymin, zmax, xmed);
    sort<z,!upz,!upx, upy> (m7, m8, zmed, xmax, ymin, zmin, xmed, ymed);
  }

  template <class RandomAccessIterator>
  void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    K k;
    double xmin=to_double(k.compute_x_3_object()(*begin)),
      ymin=to_double(k.compute_y_3_object()(*begin)),
      zmin=to_double(k.compute_z_3_object()(*begin)),
      xmax=xmin,
      ymax=ymin,
      zmax=zmin;
    for(RandomAccessIterator it=begin+1; it<end; ++it){
      if ( to_double(k.compute_x_3_object()(*it)) < xmin) 
        xmin = to_double(k.compute_x_3_object()(*it));
      if ( to_double(k.compute_y_3_object()(*it)) < ymin) 
        ymin = to_double(k.compute_y_3_object()(*it));
      if ( to_double(k.compute_z_3_object()(*it)) < zmin) 
        zmin = to_double(k.compute_z_3_object()(*it));
      if ( to_double(k.compute_x_3_object()(*it)) > xmax) 
        xmax = to_double(k.compute_x_3_object()(*it));
      if ( to_double(k.compute_y_3_object()(*it)) > ymax) 
        ymax = to_double(k.compute_y_3_object()(*it));
      if ( to_double(k.compute_z_3_object()(*it)) > zmax) 
        zmax = to_double(k.compute_z_3_object()(*it));
    }

    sort <0, false, false, false> (begin, end, xmin,ymin,zmin,xmax,ymax,zmax);
  }
};


struct Middle {};
struct Median {};


// A policy to select the sorting strategy.

template < typename Tag >
struct Hilbert_policy {};

typedef Hilbert_policy<Middle>      Hilbert_sort_middle_policy;
typedef Hilbert_policy<Median>      Hilbert_sort_median_policy;

template <class K,  class Hilbert_policy >
class Hilbert_sort_3;

template <class K>  
class Hilbert_sort_3<K, Hilbert_sort_median_policy >
  : public Hilbert_sort_median_3<K>
{
public:
  Hilbert_sort_3 (const K &k=K() , std::ptrdiff_t limit=1 )
    : Hilbert_sort_median_3<K> (k,limit)
  {}
};

template <class K>
class Hilbert_sort_3<K, Hilbert_sort_middle_policy >
  : public Hilbert_sort_middle_3<K>
{
public:
  Hilbert_sort_3 (const K &k=K() , std::ptrdiff_t limit=1 )
    : Hilbert_sort_middle_3<K> (k,limit)
  {}
};


template <class Sort>
class Multiscale_sort
{
  Sort _sort;
  std::ptrdiff_t _threshold;
  double _ratio;

public:
  Multiscale_sort (const Sort &sort = Sort(), std::ptrdiff_t threshold = 1, double ratio = 0.5)
    : _sort (sort), _threshold (threshold), _ratio (ratio)
  {
    assert (0. <= ratio && ratio <= 1.);
  }

  template <class RandomAccessIterator>
  void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    typedef typename std::iterator_traits<RandomAccessIterator>::difference_type difference_type;
    RandomAccessIterator middle = begin;
    if (end - begin >= _threshold) {
      middle = begin + difference_type ((end - begin) * _ratio);
      this->operator() (begin, middle);
    }
    _sort (middle, end);
  }
};

namespace internal {
  template <class RandomAccessIterator, class Policy, class Kernel>
  void spatial_sort (
    RandomAccessIterator begin, RandomAccessIterator end,
    const Kernel &k, 
    Policy /*policy*/,
    typename Kernel::Point_3 *,
    std::ptrdiff_t threshold_hilbert,
    std::ptrdiff_t threshold_multiscale,
    double ratio)
  {
    size_t diff=std::abs(end-begin);
    std::srand((unsigned)diff);
    typedef Hilbert_sort_3<Kernel, Policy> Sort;
    std::random_shuffle(begin,end);

    if (threshold_hilbert==0) threshold_hilbert=8;
    if (threshold_multiscale==0) threshold_multiscale=64;
    if (ratio==0.0) ratio=0.125;

    (Multiscale_sort<Sort> (Sort (k, threshold_hilbert), 
      threshold_multiscale, ratio)) (begin, end);

  }
}

template <class RandomAccessIterator, class Policy, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
  const Kernel &k,
  Policy policy,
  std::ptrdiff_t threshold_hilbert=0,
  std::ptrdiff_t threshold_multiscale=0,
  double ratio=0.0)
{
  typedef std::iterator_traits<RandomAccessIterator> ITraits;
  typedef typename ITraits::value_type               value_type;

  internal::spatial_sort(begin, end, k, policy, static_cast<value_type *> (0),
    threshold_hilbert,threshold_multiscale,ratio);
}


template <class RandomAccessIterator, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
  const Kernel &k,
  std::ptrdiff_t threshold_hilbert=0,
  std::ptrdiff_t threshold_multiscale=0,
  double ratio=0.0)
{
  spatial_sort (begin, end, k,
    Hilbert_sort_median_policy(),
    threshold_hilbert,threshold_multiscale,ratio);
}

template<class NODE_ARRAY>
class SpatialSortingDefault
{
public:
  SpatialSortingDefault(const NODE_ARRAY& nodes)
    :_nodes(nodes)
  {}

  struct Less_x_3
  {
    const NODE_ARRAY& _nodes;
    Less_x_3(const NODE_ARRAY& nodes)
      :_nodes(nodes)
    {}
    bool operator()(int p, int q) const
    {
      return _nodes[p][0] < _nodes[q][0];
    }

  };

  struct Less_y_3
  {
    const NODE_ARRAY& _nodes;
    Less_y_3(const NODE_ARRAY& nodes)
      :_nodes(nodes)
    {}
    bool operator()(int p, int q) const
    {
      return  _nodes[p][1] < _nodes[q][1];
    }
  };

  struct Less_z_3
  {
    const NODE_ARRAY& _nodes;
    Less_z_3(const NODE_ARRAY& nodes)
      :_nodes(nodes)
    {}

    bool operator()(int p, int q) const
    {
      return  _nodes[p][2] < _nodes[q][2];
    }

  };


  typedef int Point_3;

  Less_x_3 less_x_3_object() const{return Less_x_3(_nodes);}

  Less_y_3 less_y_3_object() const {return Less_y_3(_nodes);}

  Less_z_3 less_z_3_object() const  {return Less_z_3(_nodes);}

private:
  const NODE_ARRAY& _nodes;
};


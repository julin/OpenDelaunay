#pragma once
#ifdef __cplusplus
# define DELAU_NS_BEGIN namespace DELAU_NS {
# define DELAU_NS_END }
#else
#if (__BORLANDC__ <= 0x460)
typedef enum { false = 0, true } bool;
#endif
#define DELAU_NS_BEGIN
#define DELAU_NS_END
#endif /*__cplusplus*/

#include <vector>

DELAU_NS_BEGIN
void SortHilbert(const std::vector<double>&, std::vector<size_t>&);
void DelaunayTriangulation(const std::vector<double>&,std::vector<size_t> &);
DELAU_NS_END

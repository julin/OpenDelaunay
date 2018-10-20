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

#include <functional>
#include <vector>

DELAU_NS_BEGIN
typedef std::function<void(size_t, double*)> PointFunctor;

DELAU_NS_END

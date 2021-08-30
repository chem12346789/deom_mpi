#include "deom.hpp"

inline int INDEX2(const int i, const int j, const int ld) { return i * ld + j; }
inline int INDEX3(const int i, const int j, const int k, const int ld) { return i * ld * ld + j * ld + k; }
inline int INDEX4(const int i, const int j, const int k, const int l, const int ld1, const int ld2) {
  return i * ld1 * ld2 * ld2 + j * ld2 * ld2 + k * ld2 + l;
}
inline int INDEX2(const int i, const int j) { return i * NSYS + j; }
inline int INDEX3(const int i, const int j, const int k) {
  return i * NSYS * NSYS + j * NSYS + k;
}
inline int INDEX4(const int i, const int j, const int k, const int l) {
  return i * NSYS * NSYS * NSYS + j * NSYS * NSYS + k * NSYS + l;
}

inline ullint INDEX2(const ullint i, const ullint j, const ullint ld) { return i * ld + j; }
inline ullint INDEX3(const ullint i, const ullint j, const ullint k, const ullint ld) { return i * ld * ld + j * ld + k; }
inline ullint INDEX4(const ullint i, const ullint j, const ullint k, const ullint l, const ullint ld1, const ullint ld2) {
  return i * ld1 * ld2 * ld2 + j * ld2 * ld2 + k * ld2 + l;
}
inline ullint INDEX2(const ullint i, const ullint j) { return i * NSYS + j; }
inline ullint INDEX3(const ullint i, const ullint j, const ullint k) {
  return i * NSYS * NSYS + j * NSYS + k;
}
inline ullint INDEX4(const ullint i, const ullint j, const ullint k, const ullint l) {
  return i * NSYS * NSYS * NSYS + j * NSYS * NSYS + k * NSYS + l;
}

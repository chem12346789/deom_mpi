#ifndef PULSE_H_
#define PULSE_H_

#include "deom.hpp"
#include <cmath>

typedef struct {
  bool on;
  double ampl;
  double sigm;
  double freq;
  double varp;
} PULSE;

double pulse_et(double t, PULSE p) {
  double _et = 0;
  if (p.on) {
    _et = (p.ampl * 3.14159265358979323846 / (sqrt(2.0 * 3.14159265358979323846) * p.sigm)) * exp(-0.5 * t * t / (p.sigm * p.sigm)) * cos(p.freq * t + p.varp * t * t);
  }
  return _et;
}

#endif

#include "deom.hpp"

void gamma_statistics(DEOMAUX *daux, DEOM *d, const char *fname_gamma) {
  ofstream fstream_gamma(fname_gamma);
  VectorXd key(d->nind);
  for (size_t i = 0; i < d->nddo; i++) {
    key = (d->keys.row(iado)).cast<double>();
    complex<double> gams = key.transpose() * d->expn;
    fname_gamma << d->ddos[i].format(HeavyFmt) << '\n';
  }
  fname_gamma.close();
}
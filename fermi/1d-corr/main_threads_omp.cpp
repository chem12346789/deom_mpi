#include "1d_corr_common.cpp"
// #include "1d_corr_fast_slow.cpp"
#include "equilibrium.cpp"
#include "init.cpp"
#include "read_write.cpp"

void init_deom(const Json &json, int argc, char *argv[]) {
  DEOM *d, deom;
  DEOMAUX *daux, deom_aux;
  CTRL *c, ctrl;
  d = &deom;
  daux = &deom_aux;
  c = &ctrl;

  printf("init\n");
  Init_sys(d, json, "input.json");
  Init_aux(daux, d, json);
  Init_dip(d->dipole, d, json["dipole"]);
  Init_dip(d->dipole1, d, json["dipole1"]);

  Init_ctrl(c, json["equilibrium"]);
  equilibrium(daux, d, c, "prop-rho-eq1.dat");
  write_ddos(daux, d, "equ_ddos1", "equ_keys1");

  Init_ctrl(c, json);
  // read_ddos(daux, d, "equ_ddos1", "equ_keys1");
  od_corr(daux, d, c, "prop-rho1.dat", "prop-pol1.dat", 'l');
  // stochastic(daux, d, c, "prop-rho.dat");

  // read_ddos(daux, d, "tmp_ddos1", "tmp_keys1");
}
//

int main(int argc, char *argv[]) {
  ifstream jsonFile("input.json");
  stringstream strStream;
  strStream << jsonFile.rdbuf();
  string jsonStr = strStream.str();
  string err;
  const Json json = Json::parse(jsonStr, err);
  if (!err.empty()) {
    printf("Error in parsing input file: %s\n", err.c_str());
    return 1;
  }

  init_deom(json, argc, argv);
}

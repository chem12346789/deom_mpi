#include "1d_corr_fast_slow.cpp"
#include "equilibrium.cpp"
#include "init.cpp"
#include "read_write.cpp"

void init_deom(const Json &json, int argc, char *argv[]) {
  printf("init\n");
  DEOM *d, deom;
  DEOMAUX *daux, deom_aux;
  d = &deom;
  daux = &deom_aux;
  CTRL *c, ctrl;
  c = &ctrl;

  Init_sys(d, json);
  Init_aux(daux, d, json);

  Init_dip(d, json);

  FILE *frho = fopen("prop-rho.dat", "w");
  FILE *frho_eq = fopen("prop-rho-eq.dat", "w");
  FILE *fpolar = fopen("prop-pol.dat", "w");

  Init_ctrl(c, 0.05, 0, 5000);
  equilibrium(daux, d, c, frho_eq);

  write_ddos(daux, d, "tmp_ddos", "tmp_keys");
  Init_ctrl(c, json);
  od_corr(daux, d, c, frho, fpolar, 'l');
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

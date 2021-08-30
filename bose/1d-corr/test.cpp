#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;

int add(std::vector<int> &__restrict__ a, std::vector<int> &__restrict__ b) {
  a[0] = 10;
  b[0] = 12;
  return a[0] + b[0];
}
int main() {
  std::vector<int> a;
  std::vector<int> b;
  add(a, b);
}
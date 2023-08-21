
#include "matrix_new.h"
#include <helib/helib.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <helib/helib.h>
#include <vector>
#include <fstream>
#include <omp.h>

using namespace std;
using namespace NTL;

void plain_text_multiplication(vector<ZZ> M1, vector<ZZ> M2, int d)
{
   vector<vector<ZZ>> A(d,vector<ZZ>(d)), B(d,vector<ZZ>(d));
   vector<vector<ZZ>> M(d,vector<ZZ>(d));
   for (int L = 0;L<d*d;L++)
   {
    int i = L/d;
    int j = L % d;
    A[i][j] = M1[L];
    B[i][j] = M2[L];   
   }
  long unsigned int p = INT32_MAX;

   for (int i = 0 ;i<d;i++)
   {
    for (int j = 0;j<d;j++)
    {
      M[i][j] = 0;
      for (int k = 0;k<d;k++)
      {
        M[i][j] = M[i][j] + A[i][k] * B[k][j];
      }
     int a;
      if (M[i][j] % p < p/2)
        a = M[i][j] % p;
      else
        a = M[i][j] % p - p;
      std::cout << setw(8) << a << " ";
    }
    cout << endl;
   }
}

void encode_base_gen_short_int(std::vector<NTL::ZZX> &uv, helib::zzX &vdw, helib::zzX &wdud)
{
  // loading bases u_iv_j
  const int n = 5760;
  ifstream fin;
  fin.open("../../bases/Size_short_int.txt", ios::in);
  vector<long> size;
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    size.push_back(a);
  }
  fin.close();
  fin.open("../../bases/uv_coeffs_short_int.txt", ios::in);
  vector<long> uv_coeffs;
  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < size[i]; j++)
    {
      int a;
      fin >> a;
      uv_coeffs.push_back(a);
    }
  }
  fin.close();
  fin.open("../../bases/uv_deg_short_int.txt", ios::in);
  vector<long> uv_degrees;
  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < size[i]; j++)
    {
      int a;
      fin >> a;
      uv_degrees.push_back(a);
    }
  }
  fin.close();
  int l = 0;
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      uv[16 * i + j].SetLength(n);
      for (int k = 0; k < size[16 * i + j]; k++)
      {
        uv[16 * i + j][uv_degrees[l]] = uv_coeffs[l];
        l = l + 1;
      }
    }
  }
  vector<long> vdw_coeffs, wdud_coeffs;
  vector<long> vdw_degree, wdud_degree;
  ifstream fin1, fin2;
  vdw.SetLength(n);
  wdud.SetLength(n);
  for (int i = 0; i < n; i++)
  {
    vdw[i] = 0;
    wdud[i] = 0;
  }

  // loading vdw;
  int vdw_size = 614;
  fin1.open("../../bases/vdw_coeffs_short_int.txt", ios::in);
  fin2.open("../../bases/vdw_degree_short_int.txt", ios::in);
  const long inverse_17 = 30841;

  for (int i = 0; i < vdw_size; i++)
  {
    int a;
    fin1 >> a;
    vdw_coeffs.push_back(a);
    fin2 >> a;
    vdw_degree.push_back(a);
    vdw[vdw_degree[i]] = vdw_coeffs[i] * inverse_17;
  }
  fin1.close();
  fin2.close();

  // loading bases wdud;
  int wdud_size = 2510;
  fin1.open("../../bases/wdud_coeffs_short_int.txt", ios::in);
  fin2.open("../../bases/wdud_degree_short_int.txt", ios::in);
  const long inverse_475 = 55465;
  for (int i = 0; i < wdud_size; i++)
  {
    int a;
    fin1 >> a;
    wdud_coeffs.push_back(a);
    fin2 >> a;
    wdud_degree.push_back(a);
    wdud[wdud_degree[i]] = wdud_coeffs[i] * inverse_475;
  }
  fin1.close();
  fin2.close();
}

void encode_16(NTL::Vec<long> &poly, std::vector<NTL::ZZ> &M, const std::vector<NTL::ZZX> &uv)
{
  NTL::ZZX p;
  long i, j;
  for (i = 0; i < 256; i++)
  {
    p = p + M[i] * uv[i];
  }
  poly.SetLength(deg(p) + 1);
  for (i = 0; i < deg(p) + 1; i++)
  {
    NTL::conv(poly[i], p[i]);
  }
}

void decode_base_gen_short_int(std::vector<NTL::ZZX> &Base)
{
  // loading decoding bases
  const int n = 5760;
  ifstream fin;
  vector<int> M_size;
  fin.open("../../bases/M_size_short_int.txt", ios::in);
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    M_size.push_back(a);
  }
  fin.close();
  ifstream fin1, fin2;
  fin1.open("../../bases/decoding_degree_short_int.txt", ios::in);
  fin2.open("../../bases/decoding_coeffs_short_int.txt", ios::in);
  vector<int> degree;
  for (int i = 0; i < 256; i++)
  {
    int n0;
    for (int j = 0; j < M_size[i]; j++)
    {
      int s;
      fin1 >> s;
      int c;
      fin2 >> c;
      if (j == 0)
      {
        Base[i].SetLength(n);
        n0 = s;
      }
      Base[i][s] = c;
    }
    Base[i].SetLength(n0 + 1);
  }
  fin1.close();
  fin2.close();
  bool t = true;
}

void decode_16_short_int(std::vector<NTL::ZZ> &result_vector, std::vector<NTL::ZZX> &Base, NTL::ZZX f, const std::vector<NTL::ZZX> &uv)
{
  const int n = 5760;
  if (deg(f) == -1)
    return;
  vector<vector<ZZ>> decode_vectors = vector<vector<ZZ>>(256, vector<ZZ>(256));
  vector<int> combination_size(256);
  ifstream fin;
  fin.open("../../bases/combination_size_short_int.txt", ios::in);
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    combination_size[i] = a;
  }
  ifstream fin1, fin2;
  fin1.open("../../bases/decoding_combination_short_int.txt", ios::in);
  fin2.open("../../bases/combination_coeffs_short_int.txt", ios::in);

  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < combination_size[i]; j++)
    {
      int s;
      fin1 >> s;
      int c;
      fin2 >> c;
      decode_vectors[i][s - 1] = c;
    }
  }
  fin1.close();
  fin2.close();
  result_vector.resize(256);
  f.SetLength(n);
  for (int i = 0; i < 256; i++)
  {
    NTL::ZZ coeff;
    coeff = f[deg(Base[i])];
    f = f - coeff * Base[i][deg(Base[i])] * Base[i];
    for (int j = 0; j < 256; j++)
    {
      decode_vectors[i][j] =
          decode_vectors[i][j] * coeff * Base[i][deg(Base[i])];
      result_vector[j] = result_vector[j] + decode_vectors[i][j];
    }
  }
}

void HEmul_16_short_int(helib::Ctxt &result, helib::Ctxt ctxt1, helib::Ctxt ctxt2, helib::zzX &vdw, helib::zzX &wdud)
{
  // HELIB_NTIMER_START(Trace17);
  //-------------------------------
  helib::Ctxt ctxt_temp = ctxt1;
  helib::Ctxt ctxt_copy = ctxt1;
  const vector<long> Trace17 = {1901, 4276, 2376, 951};
  const vector<long> Trace19 = {851, 5526, 5101, 6801, 2551};
  const vector<long> Trace25 = {647, 6784, 7107, 324, 3231};

  // HELIB_NTIMER_START(Trace17);
  ctxt2.multByConstant(wdud);
  for (int i = 0; i < Trace17.size(); i++)
  {
    ctxt_temp = ctxt2;
    ctxt_temp.smartAutomorph(Trace17[i]);
    ctxt2.addCtxt(ctxt_temp);
  }

  //------------------------------
  // HELIB_NTIMER_STOP(Trace17);
  // helib::printNamedTimer(std::cout, "Trace17");
  //----------------------------
  // HELIB_NTIMER_START(Trace19);
  ctxt1.multByConstant(vdw);
  ctxt_copy = ctxt1;

  for (int i = 0; i < 4; i++)
  {
    ctxt_temp = ctxt1;
    ctxt_temp.smartAutomorph(Trace19[i]);
    ctxt1.addCtxt(ctxt_temp);
    if (i == 2)
    {
      ctxt_copy.smartAutomorph(Trace19[4]);
      ctxt1.addCtxt(ctxt_copy);
    }
  }
  // HELIB_NTIMER_STOP(Trace19);
  // helib::printNamedTimer(std::cout, "Trace19");
  //   //-------------------------------
  //  HELIB_NTIMER_START(mult);
  ctxt1.multiplyBy(ctxt2);
  // HELIB_NTIMER_STOP(mult);
  // helib::printNamedTimer(std::cout, "mult");
  ctxt1.dropSmallAndSpecialPrimes();
  // //-------------------------------
  // HELIB_NTIMER_START(Trace5);
  ctxt_copy = ctxt1;
  for (int i = 0; i < 4; i++)
  {
    ctxt_temp = ctxt1;
    ctxt_temp.smartAutomorph(Trace25[i]);
    ctxt1.addCtxt(ctxt_temp);
    if (i == 1)
    {
      ctxt_copy.smartAutomorph(Trace25[4]);
      ctxt1.addCtxt(ctxt_copy);
    }
  }
  // HELIB_NTIMER_STOP(Trace5);
  // helib::printNamedTimer(std::cout, "Trace5");
  result = ctxt1;
}

void HEmul_block_short_int(vector<vector<helib::Ctxt>> &ctxt, vector<vector<helib::Ctxt>> ctxt1, vector<vector<helib::Ctxt>> ctxt2, helib::zzX &vdw, helib::zzX &wdud, int M)
{
  const vector<long> Trace17 = {1901, 4276, 2376, 951};
  const vector<long> Trace19 = {851, 5526, 5101, 6801, 2551};
  const vector<long> Trace25 = {647, 6784, 7107, 324, 3231};
  // HELIB_NTIMER_START(Trace17);
  //-------------------------------

  // #pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt2[i][j];
    helib::Ctxt ctxt_copy = ctxt2[i][j];
    ctxt2[i][j].multByConstant(wdud);
    ctxt_copy = ctxt2[i][j];
    for (int k = 0; k < Trace17.size(); k++)
    {
      ctxt_temp = ctxt2[i][j];
      ctxt_temp.smartAutomorph(Trace17[k]);
      ctxt2[i][j].addCtxt(ctxt_temp);
    }
  }
  //------------------------------
  // HELIB_NTIMER_STOP(Trace27);
  // helib::printNamedTimer(std::cout, "Trace27");
  //----------------------------
  // HELIB_NTIMER_START(Trace19);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt1[i][j];
    helib::Ctxt ctxt_copy = ctxt1[i][j];
    ctxt1[i][j].multByConstant(vdw);
    ctxt_copy = ctxt1[i][j];
    for (int k = 0; k < Trace19.size() - 1; k++)
    {
      ctxt_temp = ctxt1[i][j];
      ctxt_temp.smartAutomorph(Trace19[k]);
      ctxt1[i][j].addCtxt(ctxt_temp);
      if (k == 2)
      {
        ctxt_copy.smartAutomorph(Trace19[4]);
        ctxt1[i][j].addCtxt(ctxt_copy);
      }
    }
  }
// HELIB_NTIMER_STOP(Trace19);
// helib::printNamedTimer(std::cout, "Trace17");
//-------------------------------
// HELIB_NTIMER_START(mult);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    ctxt[i][j] = ctxt1[i][0];
    ctxt[i][j].multiplyBy(ctxt2[0][j]);
    for (int k = 1; k < M; k++)
    {
      helib::Ctxt temp = ctxt1[i][k];
      temp.multiplyBy(ctxt2[k][j]);
      ctxt[i][j].addCtxt(temp);
    }
    ctxt[i][j].dropSmallAndSpecialPrimes();
  }
  // HELIB_NTIMER_STOP(mult);
  // helib::printNamedTimer(std::cout, "mult");
  //  //-------------------------------
  // HELIB_NTIMER_START(Trace25);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt[i][j];
    helib::Ctxt ctxt_copy = ctxt[i][j];
    for (int k = 0; k < Trace25.size() - 1; k++)
    {
      ctxt_temp = ctxt[i][j];
      ctxt_temp.smartAutomorph(Trace25[k]);
      ctxt[i][j].addCtxt(ctxt_temp);
      if (k == 1)
      {
        ctxt_copy.smartAutomorph(Trace25[4]);
        ctxt[i][j].addCtxt(ctxt_copy);
      }
    }
  }
  // HELIB_NTIMER_STOP(Trace25);
  // helib::printNamedTimer(std::cout, "Trace25");
}

void HEmul_16_int(helib::Ctxt &result, helib::Ctxt ctxt1, helib::Ctxt ctxt2, helib::zzX &vdw, helib::zzX &wdud)
{
  const vector<long> Trace27 = {6176, 1426, 7126, 4751, 5701};
  const vector<long> Trace19 = {8776, 3376, 4051, 10126, 7426};
  const vector<long> Trace25 = {7183, 514, 7696, 10774, 2566, 2053};
  // HELIB_NTIMER_START(Trace27);
  //-------------------------------
  helib::Ctxt ctxt_temp = ctxt1;
  helib::Ctxt ctxt_copy = ctxt1;
  ctxt2.multByConstant(wdud);
  ctxt_copy = ctxt2;
  for (int i = 0; i < Trace27.size() - 1; i++)
  {
    ctxt_temp = ctxt2;
    ctxt_temp.smartAutomorph(Trace27[i]);
    ctxt2.addCtxt(ctxt_temp);
    if (i == 2)
    {
      ctxt_temp = ctxt_copy;
      ctxt_temp.smartAutomorph(Trace27[4]);
      ctxt2.addCtxt(ctxt_temp);
    }
  }
  //------------------------------
  // HELIB_NTIMER_STOP(Trace27);
  // helib::printNamedTimer(std::cout, "Trace27");
  //----------------------------
  // HELIB_NTIMER_START(Trace19);

  ctxt1.multByConstant(vdw);
  ctxt_copy = ctxt1;
  for (int i = 0; i < Trace19.size() - 1; i++)
  {
    ctxt_temp = ctxt1;
    ctxt_temp.smartAutomorph(Trace19[i]);
    ctxt1.addCtxt(ctxt_temp);
    if (i == 2)
    {
      ctxt_temp = ctxt_copy;
      ctxt_temp.smartAutomorph(Trace19[4]);
      ctxt1.addCtxt(ctxt_temp);
    }
  }

  // HELIB_NTIMER_STOP(Trace19);
  // helib::printNamedTimer(std::cout, "Trace19");
  //   //-------------------------------
  // HELIB_NTIMER_START(mult);
  ctxt1.multiplyBy(ctxt2);
  // HELIB_NTIMER_STOP(mult);
  // helib::printNamedTimer(std::cout, "mult");
  ctxt1.dropSmallAndSpecialPrimes();
  // //-------------------------------
  // HELIB_NTIMER_START(Trace5);
  ctxt_copy = ctxt1;
  for (int i = 0; i < Trace25.size() - 2; i++)
  {
    ctxt_temp = ctxt1;
    ctxt_temp.smartAutomorph(Trace25[i]);
    ctxt1.addCtxt(ctxt_temp);
    if (i == 2)
    {
      ctxt_temp = ctxt_copy;
      ctxt_temp.smartAutomorph(Trace25[4]);
      ctxt1.addCtxt(ctxt_temp);
      ctxt_temp = ctxt_copy;
      ctxt_temp.smartAutomorph(Trace25[5]);
      ctxt1.addCtxt(ctxt_temp);
    }
  }
  // HELIB_NTIMER_STOP(Trace5);
  // helib::printNamedTimer(std::cout, "Trace5");
  result = ctxt1;
}

void HEmul_block_int(vector<vector<helib::Ctxt>> &ctxt, vector<vector<helib::Ctxt>> ctxt1, vector<vector<helib::Ctxt>> ctxt2, helib::zzX &vdw, helib::zzX &wdud, int M)
{
  const vector<long> Trace27 = {1426, 6176, 7126, 4751, 5701};
  const vector<long> Trace19 = {8776, 3376, 4051, 10126, 7426};
  const vector<long> Trace25 = {1027, 7183, 7696, 10774, 2566, 2053};
  // HELIB_NTIMER_START(Trace27);
  //-------------------------------

#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt2[i][j];
    helib::Ctxt ctxt_copy = ctxt2[i][j];
    ctxt2[i][j].multByConstant(wdud);
    ctxt_copy = ctxt2[i][j];
    for (int k = 0; k < Trace27.size() - 1; k++)
    {
      ctxt_temp = ctxt2[i][j];
      ctxt_temp.smartAutomorph(Trace27[k]);
      ctxt2[i][j].addCtxt(ctxt_temp);
      if (k == 2)
      {
        ctxt_temp = ctxt_copy;
        ctxt_temp.smartAutomorph(Trace27[4]);
        ctxt2[i][j].addCtxt(ctxt_temp);
      }
    }
  }
//------------------------------
// HELIB_NTIMER_STOP(Trace27);
// helib::printNamedTimer(std::cout, "Trace27");
//----------------------------
// HELIB_NTIMER_START(Trace19);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt1[i][j];
    helib::Ctxt ctxt_copy = ctxt1[i][j];
    ctxt1[i][j].multByConstant(vdw);
    ctxt_copy = ctxt1[i][j];
    for (int k = 0; k < Trace19.size() - 1; k++)
    {
      ctxt_temp = ctxt1[i][j];
      ctxt_temp.smartAutomorph(Trace19[k]);
      ctxt1[i][j].addCtxt(ctxt_temp);
      if (k == 2)
      {
        ctxt_temp = ctxt_copy;
        ctxt_temp.smartAutomorph(Trace19[4]);
        ctxt1[i][j].addCtxt(ctxt_temp);
      }
    }
  }
// HELIB_NTIMER_STOP(Trace19);
// helib::printNamedTimer(std::cout, "Trace19");
//   //-------------------------------
// HELIB_NTIMER_START(mult);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    ctxt[i][j] = ctxt1[i][0];
    ctxt[i][j].multiplyBy(ctxt2[0][j]);
    for (int k = 1; k < M; k++)
    {
      helib::Ctxt temp = ctxt1[i][k];
      temp.multiplyBy(ctxt2[k][j]);
      ctxt[i][j].addCtxt(temp);
    }
    ctxt[i][j].dropSmallAndSpecialPrimes();
  }
  // HELIB_NTIMER_STOP(mult);
  // helib::printNamedTimer(std::cout, "mult");
  //  //-------------------------------
  HELIB_NTIMER_START(Trace5);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt[i][j];
    helib::Ctxt ctxt_copy = ctxt[i][j];
    for (int k = 0; k < Trace25.size() - 2; k++)
    {
      ctxt_temp = ctxt[i][j];
      ctxt_temp.smartAutomorph(Trace25[k]);
      ctxt[i][j].addCtxt(ctxt_temp);
      if (k == 2)
      {
        ctxt_temp = ctxt_copy;
        ctxt_temp.smartAutomorph(Trace25[4]);
        ctxt[i][j].addCtxt(ctxt_temp);
        ctxt_temp = ctxt_copy;
        ctxt_temp.smartAutomorph(Trace25[5]);
        ctxt[i][j].addCtxt(ctxt_temp);
      }
    }
  }
  // HELIB_NTIMER_STOP(Trace5);
  // helib::printNamedTimer(std::cout, "Trace5");
}

void printf_matrix(const std::vector<NTL::ZZ> &M, int n, int w, long unsigned int &p)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int a;
      if (M[i * n + j] % p < p/2)
        a = M[i * n + j] % p;
      else
        a = M[i * n + j] % p - p;
      std::cout << setw(w) << a << " ";
    }
    std::cout << std::endl;
  }
}

void printf_matrix(const vector<ZZ> &M, int n, int w, long unsigned int &p, ofstream &fout)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      int a;
      if (M[i * n + j] % p < p/2)
        a = M[i * n + j] % p;
      else
        a = M[i * n + j] % p - p;
      fout << setw(w) << a << " ";
    }
    fout << std::endl;
  }
  cout << "   See output" << endl;
}

void encode_base_gen_16_int(std::vector<NTL::ZZX> &uv, helib::zzX &vdw, helib::zzX &wdud)
{
  // loading bases u_iv_j
  const int n = 6480;
  ifstream fin;
  fin.open("../../bases/Size.txt", ios::in);
  vector<long> size;
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    size.push_back(a);
  }
  fin.close();
  fin.open("../../bases/v_coeffs.txt", ios::in);
  vector<long> uv_coeffs;
  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < size[i]; j++)
    {
      int a;
      fin >> a;
      uv_coeffs.push_back(a);
    }
  }
  fin.close();
  fin.open("../../bases/v_deg.txt", ios::in);
  vector<long> uv_degrees;
  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < size[i]; j++)
    {
      int a;
      fin >> a;
      uv_degrees.push_back(a);
    }
  }
  fin.close();
  int l = 0;
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      uv[16 * i + j].SetLength(n);
      for (int k = 0; k < size[16 * i + j]; k++)
      {
        uv[16 * i + j][uv_degrees[l]] = uv_coeffs[l];
        l = l + 1;
      }
    }
  }

  vector<long> vdw_coeffs, wdud_coeffs;
  vector<long> vdw_degree, wdud_degree;
  ifstream fin1, fin2;
  vdw.SetLength(n);
  wdud.SetLength(n);
  for (int i = 0; i < n; i++)
  {
    vdw[i] = 0;
    wdud[i] = 0;
  }

  // loading vdw;
  int vdw_size = 85;
  fin1.open("../../bases/vdw_coeffs.txt", ios::in);
  fin2.open("../../bases/vdw_degree.txt", ios::in);
  const long inverse_19 = 1017229096;

  for (int i = 0; i < vdw_size; i++)
  {
    int a;
    fin1 >> a;
    vdw_coeffs.push_back(a);
    fin2 >> a;
    vdw_degree.push_back(a);
    vdw[vdw_degree[i]] = vdw_coeffs[i] * inverse_19;
  }
  fin1.close();
  fin2.close();

  // loading bases wdud;
  int wdud_size = 242;
  fin1.open("../../bases/wdud_coeffs.txt", ios::in);
  fin2.open("../../bases/wdud_degree.txt", ios::in);
  const long inverse_675 = 1485740538;
  for (int i = 0; i < wdud_size; i++)
  {
    int a;
    fin1 >> a;
    wdud_coeffs.push_back(a);
    fin2 >> a;
    wdud_degree.push_back(a);
    wdud[wdud_degree[i]] = wdud_coeffs[i] * inverse_675;
  }
  fin1.close();
  fin2.close();
}

void decode_base_gen_16_int(std::vector<NTL::ZZX> &Base)
{
  // loading decoding bases
  const int n = 6480;
  ifstream fin;
  vector<int> M_size;
  fin.open("../../bases/M_size.txt", ios::in);
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    M_size.push_back(a);
  }
  fin.close();
  ifstream fin1, fin2;
  fin1.open("../../bases/decoding_degree.txt", ios::in);
  fin2.open("../../bases/decoding_coeffs.txt", ios::in);
  vector<int> degree;
  for (int i = 0; i < 256; i++)
  {
    int n0;
    for (int j = 0; j < M_size[i]; j++)
    {
      int s;
      fin1 >> s;
      int c;
      fin2 >> c;
      if (j == 0)
      {
        Base[i].SetLength(n);
        n0 = s;
      }
      Base[i][s] = c;
    }
    Base[i].SetLength(n0 + 1);
  }
  fin1.close();
  fin2.close();
  bool t = true;
}

void decode_16_int(std::vector<NTL::ZZ> &result_vector, std::vector<NTL::ZZX> &Base, NTL::ZZX f, const std::vector<NTL::ZZX> &uv)
{
  // HELIB_NTIMER_START(decode123);
  const int n = 6480;
  if (deg(f) == -1)
    return;
  vector<vector<ZZ>> decode_vectors = vector<vector<ZZ>>(256, vector<ZZ>(256));
  vector<int> combination_size(256);
  ifstream fin;
  fin.open("../../bases/combination_size.txt", ios::in);
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    combination_size[i] = a;
  }
  ifstream fin1, fin2;
  fin1.open("../../bases/decoding_combination.txt", ios::in);
  fin2.open("../../bases/combination_coeffs.txt", ios::in);

  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < combination_size[i]; j++)
    {
      int s;
      fin1 >> s;
      int c;
      fin2 >> c;
      decode_vectors[i][s - 1] = c;
    }
  }
  fin1.close();
  fin2.close();
  result_vector.resize(256);
  f.SetLength(n);
  for (int i = 0; i < 256; i++)
  {
    NTL::ZZ coeff;
    coeff = f[deg(Base[i])];
    f = f - coeff * Base[i][deg(Base[i])] * Base[i];
    for (int j = 0; j < 256; j++)
    {
      decode_vectors[i][j] =
          decode_vectors[i][j] * coeff * Base[i][deg(Base[i])];
      result_vector[j] = result_vector[j] + decode_vectors[i][j];
    }
  }
  // cout << f << endl;
  // HELIB_NTIMER_STOP(decode123);
  // helib::printNamedTimer(std::cout, "decode123");
}

void encode_base_gen_16_int_depth_4(std::vector<NTL::ZZX> &uv, helib::zzX &vdw, helib::zzX &wdud)
{
  // loading bases u_iv_j
  const int n = 17280;
  ifstream fin;
  fin.open("../../bases/Size_depth_4.txt", ios::in);
  vector<long> size;
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    size.push_back(a);
  }
  fin.close();
  fin.open("../../bases/v_coeffs_depth_4.txt", ios::in);
  vector<long> uv_coeffs;
  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < size[i]; j++)
    {
      int a;
      fin >> a;
      uv_coeffs.push_back(a);
    }
  }
  fin.close();
  fin.open("../../bases/v_deg_depth_4.txt", ios::in);
  vector<long> uv_degrees;
  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < size[i]; j++)
    {
      int a;
      fin >> a;
      uv_degrees.push_back(a);
    }
  }
  fin.close();
  int l = 0;
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      uv[16 * i + j].SetLength(n);
      for (int k = 0; k < size[16 * i + j]; k++)
      {
        uv[16 * i + j][uv_degrees[l]] = uv_coeffs[l];
        l = l + 1;
      }
    }
  }

  vector<long> vdw_coeffs, wdud_coeffs;
  vector<long> vdw_degree, wdud_degree;
  ifstream fin1, fin2;
  vdw.SetLength(n);
  wdud.SetLength(n);
  for (int i = 0; i < n; i++)
  {
    vdw[i] = 0;
    wdud[i] = 0;
  }

  // loading vdw;
  int vdw_size = 80;
  fin1.open("../../bases/vdw_coeffs_depth_4.txt", ios::in);
  fin2.open("../../bases/vdw_degree_depth_4.txt", ios::in);
  const long inverse_17 = 252645135;

  for (int i = 0; i < vdw_size; i++)
  {
    int a;
    fin1 >> a;
    vdw_coeffs.push_back(a);
    fin2 >> a;
    vdw_degree.push_back(a);
    vdw[vdw_degree[i]] = vdw_coeffs[i] * inverse_17;
  }
  fin1.close();
  fin2.close();

  // loading bases wdud;
  int wdud_size = 265;
  fin1.open("../../bases/wdud_coeffs_depth_4.txt", ios::in);
  fin2.open("../../bases/wdud_degree_depth_4.txt", ios::in);
  const long inverse_2025 = 495246846;
  for (int i = 0; i < wdud_size; i++)
  {
    int a;
    fin1 >> a;
    wdud_coeffs.push_back(a);
    fin2 >> a;
    wdud_degree.push_back(a);
    wdud[wdud_degree[i]] = wdud_coeffs[i] * inverse_2025;
  }
  fin1.close();
  fin2.close();
}

void decode_base_gen_16_int_depth_4(std::vector<NTL::ZZX> &Base)
{
  // loading decoding bases
  const int n = 17280;
  ifstream fin;
  vector<int> M_size;
  fin.open("../../bases/M_size_depth_4.txt", ios::in);
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    M_size.push_back(a);
  }
  fin.close();
  ifstream fin1, fin2;
  fin1.open("../../bases/decoding_degree_depth_4.txt", ios::in);
  fin2.open("../../bases/decoding_coeffs_depth_4.txt", ios::in);
  vector<int> degree;
  for (int i = 0; i < 256; i++)
  {
    int n0;
    for (int j = 0; j < M_size[i]; j++)
    {
      int s;
      fin1 >> s;
      int c;
      fin2 >> c;
      if (j == 0)
      {
        Base[i].SetLength(n);
        n0 = s;
      }
      Base[i][s] = c;
    }
    Base[i].SetLength(n0 + 1);
  }
  fin1.close();
  fin2.close();
  bool t = true;
}

void decode_16_int_depth_4(std::vector<NTL::ZZ> &result_vector, std::vector<NTL::ZZX> &Base, NTL::ZZX f, const std::vector<NTL::ZZX> &uv)
{
  // HELIB_NTIMER_START(decode123);
  const int n = 17280;
  if (deg(f) == -1)
    return;
  vector<vector<ZZ>> decode_vectors = vector<vector<ZZ>>(256, vector<ZZ>(256));
  vector<int> combination_size(256);
  ifstream fin;
  fin.open("../../bases/combination_size_depth_4.txt", ios::in);
  for (int i = 0; i < 256; i++)
  {
    int a;
    fin >> a;
    combination_size[i] = a;
  }
  ifstream fin1, fin2;
  fin1.open("../../bases/decoding_combination_depth_4.txt", ios::in);
  fin2.open("../../bases/combination_coeffs_depth_4.txt", ios::in);

  for (int i = 0; i < 256; i++)
  {
    for (int j = 0; j < combination_size[i]; j++)
    {
      int s;
      fin1 >> s;
      int c;
      fin2 >> c;
      decode_vectors[i][s - 1] = c;
    }
  }
  fin1.close();
  fin2.close();
  result_vector.resize(256);
  f.SetLength(n);
  for (int i = 0; i < 256; i++)
  {
    NTL::ZZ coeff;
    coeff = f[deg(Base[i])];
    f = f - coeff * Base[i][deg(Base[i])] * Base[i];
    for (int j = 0; j < 256; j++)
    {
      decode_vectors[i][j] =
          decode_vectors[i][j] * coeff * Base[i][deg(Base[i])];
      result_vector[j] = result_vector[j] + decode_vectors[i][j];
    }
  }
  // cout << f << endl;
  // HELIB_NTIMER_STOP(decode123);
  // helib::printNamedTimer(std::cout, "decode123");
}

void HEmul_16_int_depth_4(helib::Ctxt &result, helib::Ctxt ctxt1, helib::Ctxt ctxt2, helib::zzX &vdw, helib::zzX &wdud)
{
  // g, g^2, g^4, g^9, g^8;
  const vector<long> Trace27 = {2126, 10201, 5101, 3401, 6376};
  // g,g^2,g^4,g^8;
  const vector<long> Trace17 = {2026, 8101, 12151, 32401};
  // g,g^2,g^4,g^10,g^8,g^9;
  const vector<long> Trace25 = {1378, 5509, 20656, 33049, 6886, 22033};
  // HELIB_NTIMER_START(Trace27);
  //-------------------------------
  helib::Ctxt ctxt_temp = ctxt1;
  helib::Ctxt ctxt_copy = ctxt1;
  ctxt2.multByConstant(wdud);
  ctxt_copy = ctxt2;
  for (int i = 0; i < Trace27.size() - 1; i++)
  {
    ctxt_temp = ctxt2;
    ctxt_temp.smartAutomorph(Trace27[i]);
    ctxt2.addCtxt(ctxt_temp);
    if (i == 2)
    {
      ctxt_temp = ctxt_copy;
      ctxt_temp.smartAutomorph(Trace27[4]);
      ctxt2.addCtxt(ctxt_temp);
    }
  }
  //------------------------------
  // HELIB_NTIMER_STOP(Trace27);
  // helib::printNamedTimer(std::cout, "Trace27");
  //----------------------------
  // HELIB_NTIMER_START(Trace17);

  ctxt1.multByConstant(vdw);
  ctxt_copy = ctxt1;
  for (int i = 0; i < Trace17.size(); i++)
  {
    ctxt_temp = ctxt1;
    ctxt_temp.smartAutomorph(Trace17[i]);
    ctxt1.addCtxt(ctxt_temp);
  }

  // HELIB_NTIMER_STOP(Trace19);
  // helib::printNamedTimer(std::cout, "Trace17");
  //   //-------------------------------
  // HELIB_NTIMER_START(mult);
  ctxt1.multiplyBy(ctxt2);
  // HELIB_NTIMER_STOP(mult);
  // helib::printNamedTimer(std::cout, "mult");
  ctxt1.dropSmallAndSpecialPrimes();
  // //-------------------------------
  // HELIB_NTIMER_START(Trace5);
  ctxt_copy = ctxt1;
  for (int i = 0; i < Trace25.size() - 2; i++)
  {
    ctxt_temp = ctxt1;
    ctxt_temp.smartAutomorph(Trace25[i]);
    ctxt1.addCtxt(ctxt_temp);
    if (i == 2)
    {
      ctxt_temp = ctxt_copy;
      ctxt_temp.smartAutomorph(Trace25[4]);
      ctxt1.addCtxt(ctxt_temp);
      ctxt_temp = ctxt_copy;
      ctxt_temp.smartAutomorph(Trace25[5]);
      ctxt1.addCtxt(ctxt_temp);
    }
  }
  // HELIB_NTIMER_STOP(Trace5);
  // helib::printNamedTimer(std::cout, "Trace5");
  result = ctxt1;
}

void HEmul_block_int_depth_4(vector<vector<helib::Ctxt>> &ctxt, vector<vector<helib::Ctxt>> ctxt1, vector<vector<helib::Ctxt>> ctxt2, helib::zzX &vdw, helib::zzX &wdud, int M)
{
  // g, g^2, g^4, g^9, g^8;
  const vector<long> Trace27 = {2126, 10201, 5101, 3401, 6376};
  // g,g^2,g^4,g^8;
  const vector<long> Trace17 = {2026, 8101, 12151, 32401};
  // g,g^2,g^4,g^10,g^8,g^9;
  const vector<long> Trace25 = {1378, 5509, 20656, 33049, 6886, 22033};
  // HELIB_NTIMER_START(Trace27);
  //-------------------------------

#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt2[i][j];
    helib::Ctxt ctxt_copy = ctxt2[i][j];
    ctxt2[i][j].multByConstant(wdud);
    ctxt_copy = ctxt2[i][j];
    for (int k = 0; k < Trace27.size() - 1; k++)
    {
      ctxt_temp = ctxt2[i][j];
      ctxt_temp.smartAutomorph(Trace27[k]);
      ctxt2[i][j].addCtxt(ctxt_temp);
      if (k == 2)
      {
        ctxt_temp = ctxt_copy;
        ctxt_temp.smartAutomorph(Trace27[4]);
        ctxt2[i][j].addCtxt(ctxt_temp);
      }
    }
  }
//------------------------------
// HELIB_NTIMER_STOP(Trace27);
// helib::printNamedTimer(std::cout, "Trace27");
//----------------------------
// HELIB_NTIMER_START(Trace17);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt1[i][j];
    helib::Ctxt ctxt_copy = ctxt1[i][j];
    ctxt1[i][j].multByConstant(vdw);
    ctxt_copy = ctxt1[i][j];
    for (int k = 0; k < Trace17.size(); k++)
    {
      ctxt_temp = ctxt1[i][j];
      ctxt_temp.smartAutomorph(Trace17[k]);
      ctxt1[i][j].addCtxt(ctxt_temp);
    }
  }
// HELIB_NTIMER_STOP(Trace17);
// helib::printNamedTimer(std::cout, "Trace17");
//   //-------------------------------
// HELIB_NTIMER_START(mult);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    ctxt[i][j] = ctxt1[i][0];
    ctxt[i][j].multiplyBy(ctxt2[0][j]);
    for (int k = 1; k < M; k++)
    {
      helib::Ctxt temp = ctxt1[i][k];
      temp.multiplyBy(ctxt2[k][j]);
      ctxt[i][j].addCtxt(temp);
    }
    ctxt[i][j].dropSmallAndSpecialPrimes();
  }
  // HELIB_NTIMER_STOP(mult);
  // helib::printNamedTimer(std::cout, "mult");
  //  //-------------------------------
  // HELIB_NTIMER_START(Trace25);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    helib::Ctxt ctxt_temp = ctxt[i][j];
    helib::Ctxt ctxt_copy = ctxt[i][j];
    for (int k = 0; k < Trace25.size() - 2; k++)
    {
      ctxt_temp = ctxt[i][j];
      ctxt_temp.smartAutomorph(Trace25[k]);
      ctxt[i][j].addCtxt(ctxt_temp);
      if (k == 2)
      {
        ctxt_temp = ctxt_copy;
        ctxt_temp.smartAutomorph(Trace25[4]);
        ctxt[i][j].addCtxt(ctxt_temp);
        ctxt_temp = ctxt_copy;
        ctxt_temp.smartAutomorph(Trace25[5]);
        ctxt[i][j].addCtxt(ctxt_temp);
      }
    }
  }
  // HELIB_NTIMER_STOP(Trace25);
  // helib::printNamedTimer(std::cout, "Trace25");
}

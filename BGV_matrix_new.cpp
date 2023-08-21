#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <helib/helib.h>
#include "matrix_new.h"
#include <vector>
#include <omp.h>

using namespace std;
using namespace NTL;

void matrix_16_short_int()
{
  int m1 = 17;
  int m2 = 19;
  int m3 = 25;
  unsigned long p = 65537;
  unsigned long m = m1 * m2 * m3;
  unsigned long r = 1;
  unsigned long bits = 110;
  unsigned long c = 2;

  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m1 << "*" << m2 << "*" << m3
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);

  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace17 = {1901, 4276, 2376, 951};
  const vector<long> Trace19 = {851, 5526, 5101, 6801, 2551};
  const vector<long> Trace25 = {647, 6784, 7107, 324, 3231};
  for (int i = 0; i < Trace17.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace17[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_short_int(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_short_int(Base);
  // Read matrices M1 and M2;
  vector<ZZ> M1(256), M2(256);
  ifstream fin1, fin2;
  fin1.open("../../input/A_16.txt", ios::in);
  fin2.open("../../input/B_16.txt", ios::in);
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      fin1 >> M1[i * 16 + j];
      fin2 >> M2[i * 16 + j];
    }
  }
  fin1.close();
  fin2.close();
  cout << endl;
  cout << "The matrix A is equal to: " << endl;
  printf_matrix(M1, 16, 7, p);
  cout << endl;
  cout << "The matrix B is equal to: " << endl;
  printf_matrix(M2, 16, 7, p);
  cout << endl;
  // cout << "Plaintext multiplication: A*B is equal to: " << endl;
  // plain_text_multiplication(M1,M2,16);
  // cout << endl;

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  NTL::Vec<long> poly1, poly2;
  HELIB_NTIMER_START(Encode);
  encode_16(poly1, M1, uv);
  encode_16(poly2, M2, uv);

  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  helib::Ctxt ctxt1(public_key);
  public_key.Encrypt(ctxt1, poly1);
  helib::Ctxt ctxt2(public_key);
  public_key.Encrypt(ctxt2, poly2);
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");

  // Matrix addition:
  helib::Ctxt ctxt_add = ctxt1;
  HELIB_NTIMER_START(MatAdd);
  ctxt_add.addCtxt(ctxt2);
  HELIB_NTIMER_STOP(MatAdd);
  helib::printNamedTimer(std::cout, "MatAdd");

  // Decrypt:
  helib::Ptxt<helib::BGV> plaintext_result_add(context);
  HELIB_NTIMER_START(Decrypt_add);
  secret_key.Decrypt(plaintext_result_add, ctxt_add);
  HELIB_NTIMER_STOP(Decrypt_add);
  helib::printNamedTimer(std::cout, "Decrypt_add");

  // Decode:
  HELIB_NTIMER_START(Decode_add);
  vector<ZZ> result_add(256);
  decode_16_short_int(result_add, Base, plaintext_result_add.getPolyRepr(), uv);
  HELIB_NTIMER_STOP(Decode_add);
  helib::printNamedTimer(std::cout, "Decode_add");

  cout << endl;

  cout << "Ciphertext addition result: The matrix A + B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  printf_matrix(result_add, 16, 7, p);

  cout << endl;

  cout << "The running times (in seconds) is: " << endl;
  // Matrix Multiplication:
  helib::Ctxt ctxt = ctxt1;
  HELIB_NTIMER_START(MatMult);
  HEmul_16_short_int(ctxt, ctxt1, ctxt2, vdw, wdud);
  HELIB_NTIMER_STOP(MatMult);
  helib::printNamedTimer(std::cout, "MatMult");

  // Decrypt:
  helib::Ptxt<helib::BGV> plaintext_result(context);
  HELIB_NTIMER_START(Decrypt_mult);
  secret_key.Decrypt(plaintext_result, ctxt);
  HELIB_NTIMER_STOP(Decrypt_mult);
  helib::printNamedTimer(std::cout, "Decrypt_mult");

  // Decode:
  HELIB_NTIMER_START(Decode_mult);
  vector<ZZ> result(256);
  decode_16_short_int(result, Base, plaintext_result.getPolyRepr(), uv);
  HELIB_NTIMER_STOP(Decode_mult);
  helib::printNamedTimer(std::cout, "Decode_mult");

  cout << endl;

  cout << "Ciphertext multiplication result: The matrix A * B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  printf_matrix(result, 16, 7, p);
}

void matrix_32_short_int(int thread)
{
  omp_set_num_threads(thread);

  int M = 2;
  int d = 32;

  unsigned long p = 65537;
  unsigned long m = 17 * 19 * 25;
  unsigned long r = 1;
  unsigned long bits = 110;
  unsigned long c = 2;

  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  nslots n = "
            << "5760"
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);

  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace17 = {1901, 4276, 2376, 951};
  const vector<long> Trace19 = {851, 5526, 5101, 6801, 2551};
  const vector<long> Trace25 = {647, 6784, 7107, 324, 3231};
  for (int i = 0; i < Trace17.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace17[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_short_int(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_short_int(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(d * d), M2(d * d);
  ifstream fin1, fin2;
  fin1.open("../../input/A_32.txt", ios::in);
  fin2.open("../../input/B_32.txt", ios::in);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      fin1 >> M1[i * d + j];
      fin2 >> M2[i * d + j];
    }
  }
  fin1.close();
  fin2.close();
  cout << endl;
  cout << "The matrix A is equal to: " << endl;
  printf_matrix(M1, d, 5, p);
  cout << endl;
  cout << "The matrix B is equal to: " << endl;
  printf_matrix(M2, d, 5, p);
  cout << endl;

  // cout << "Plaintext multiplication: A*B is equal to: " << endl;
  // plain_text_multiplication(M1,M2,32);
  // cout << endl;

  // Block matrices
  vector<vector<vector<ZZ>>> b_M1(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  vector<vector<vector<ZZ>>> b_M2(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  for (int i = 0; i < d * d; i++)
  {
    int i0 = i / d;
    int j0 = i % d;
    b_M1[i0 / 16][j0 / 16].push_back(M1[i]);
    b_M2[i0 / 16][j0 / 16].push_back(M2[i]);
  }

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  vector<vector<NTL::Vec<long>>> poly1(M, vector<NTL::Vec<long>>(M)), poly2(M, vector<NTL::Vec<long>>(M));
  HELIB_NTIMER_START(Encode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    encode_16(poly1[i][j], b_M1[i][j], uv);
    encode_16(poly2[i][j], b_M2[i][j], uv);
  }
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  vector<vector<helib::Ctxt>> ctxt1(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  vector<vector<helib::Ctxt>> ctxt2(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    public_key.Encrypt(ctxt1[i][j], poly1[i][j]);
    public_key.Encrypt(ctxt2[i][j], poly2[i][j]);
  }
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");

  // Matrix addition:
  vector<vector<helib::Ctxt>> ctxt_add(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HELIB_NTIMER_START(MatAdd);
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < M; j++)
    {
      ctxt_add[i][j] = ctxt1[i][j];
      ctxt_add[i][j].addCtxt(ctxt2[i][j]);
    }
  }
  HELIB_NTIMER_STOP(MatAdd);
  helib::printNamedTimer(std::cout, "MatAdd");

  HELIB_NTIMER_START(Decrypt_add);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result_add(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result_add[i][j], ctxt_add[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt_add);
  helib::printNamedTimer(std::cout, "Decrypt_add");

  // Decode:
  vector<vector<vector<ZZ>>> result_add(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode_add);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_short_int(result_add[i][j], Base, plaintext_result_add[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode_add);
  helib::printNamedTimer(std::cout, "Decode_add");

  cout << endl;

  cout << "The matrix A + B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total_add(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total_add[i * d + j] = result_add[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  ofstream fout;
  fout.open("./M_32_add.txt", ios::out);
  printf_matrix(result_total_add, d, 11, p, fout);
  cout << "   See M_32_add.txt" << endl;
  fout.close();

  cout << endl;
  cout << "The running times (in seconds) is: " << endl;
  // // Matrix Multiplication:
  HELIB_NTIMER_START(MatMult);
  vector<vector<helib::Ctxt>> ctxt(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HEmul_block_short_int(ctxt, ctxt1, ctxt2, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult);
  helib::printNamedTimer(std::cout, "MatMult");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt_mult);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result[i][j], ctxt[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt_mult);
  helib::printNamedTimer(std::cout, "Decrypt_mult");

  // Decode:
  vector<vector<vector<ZZ>>> result(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode_mult);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_short_int(result[i][j], Base, plaintext_result[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode_mult);
  helib::printNamedTimer(std::cout, "Decode_mult");

  cout << endl;

  cout << "The matrix A*B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total[i * d + j] = result[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  fout.open("./M_32_mult.txt", ios::out);
  printf_matrix(result_total, d, 11, p, fout);
  cout << "   See M_32_mult.txt" << endl;
  fout.close();
}

void matrix_64_short_int(int thread = 16)
{
  omp_set_num_threads(thread);

  int M = 4;
  int d = 64;

  unsigned long p = 65537;
  unsigned long m = 17 * 19 * 25;
  unsigned long r = 1;
  unsigned long bits = 110;
  unsigned long c = 2;

  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  nslots n = "
            << "5760"
            << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);

  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace17 = {1901, 4276, 2376, 951};
  const vector<long> Trace19 = {851, 5526, 5101, 6801, 2551};
  const vector<long> Trace25 = {647, 6784, 7107, 324, 3231};
  for (int i = 0; i < Trace17.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace17[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_short_int(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_short_int(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(d * d), M2(d * d);
  ifstream fin1, fin2;
  fin1.open("../../input/A_64.txt", ios::in);
  fin2.open("../../input/B_64.txt", ios::in);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      fin1 >> M1[i * d + j];
      fin2 >> M2[i * d + j];
    }
  }
  fin1.close();
  fin2.close();
  cout << endl;
  cout << "The matrix A is equal to: " << endl;
  printf_matrix(M1, d, 5, p);
  cout << endl;
  cout << "The matrix B is equal to: " << endl;
  printf_matrix(M2, d, 5, p);
  cout << endl;

  // cout << "Plaintext multiplication: A*B is equal to: " << endl;
  // plain_text_multiplication(M1,M2,d);
  // cout << endl;
  // Block matrices
  vector<vector<vector<ZZ>>> b_M1(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  vector<vector<vector<ZZ>>> b_M2(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  for (int i = 0; i < d * d; i++)
  {
    int i0 = i / d;
    int j0 = i % d;
    b_M1[i0 / 16][j0 / 16].push_back(M1[i]);
    b_M2[i0 / 16][j0 / 16].push_back(M2[i]);
  }

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  vector<vector<NTL::Vec<long>>> poly1(M, vector<NTL::Vec<long>>(M)), poly2(M, vector<NTL::Vec<long>>(M));
  HELIB_NTIMER_START(Encode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    encode_16(poly1[i][j], b_M1[i][j], uv);
    encode_16(poly2[i][j], b_M2[i][j], uv);
  }
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  vector<vector<helib::Ctxt>> ctxt1(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  vector<vector<helib::Ctxt>> ctxt2(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    public_key.Encrypt(ctxt1[i][j], poly1[i][j]);
    public_key.Encrypt(ctxt2[i][j], poly2[i][j]);
  }
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");

  // Matrix addition:
  vector<vector<helib::Ctxt>> ctxt_add(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HELIB_NTIMER_START(MatAdd);
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < M; j++)
    {
      ctxt_add[i][j] = ctxt1[i][j];
      ctxt_add[i][j].addCtxt(ctxt2[i][j]);
    }
  }
  HELIB_NTIMER_STOP(MatAdd);
  helib::printNamedTimer(std::cout, "MatAdd");

  HELIB_NTIMER_START(Decrypt_add);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result_add(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result_add[i][j], ctxt_add[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt_add);
  helib::printNamedTimer(std::cout, "Decrypt_add");

  // Decode:
  vector<vector<vector<ZZ>>> result_add(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode_add);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_short_int(result_add[i][j], Base, plaintext_result_add[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode_add);
  helib::printNamedTimer(std::cout, "Decode_add");

  cout << endl;

  cout << "The matrix A + B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total_add(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total_add[i * d + j] = result_add[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  ofstream fout;
  fout.open("./M_64_add.txt", ios::out);
  printf_matrix(result_total_add, d, 11, p, fout);
  cout << "   See M_64_add.txt" << endl;
  fout.close();

  cout << endl;

  cout << "The running times (in seconds) is: " << endl;

  // // Matrix Multiplication:
  HELIB_NTIMER_START(MatMult_new);
  vector<vector<helib::Ctxt>> ctxt(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HEmul_block_short_int(ctxt, ctxt1, ctxt2, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_new);
  helib::printNamedTimer(std::cout, "MatMult_new");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt_mult);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result[i][j], ctxt[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt_mult);
  helib::printNamedTimer(std::cout, "Decrypt_mult");

  // Decode:
  vector<vector<vector<ZZ>>> result(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode_mult);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_short_int(result[i][j], Base, plaintext_result[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode_mult);
  helib::printNamedTimer(std::cout, "Decode_mult");

  cout << endl;

  cout << "The matrix A*B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total[i * d + j] = result[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  fout.open("./M_64_mult.txt", ios::out);
  printf_matrix(result_total, d, 11, p, fout);
  cout << "   See M_64_mult.txt" << endl;
  fout.close();
}

void matrix_16_int()
{
  int m1 = 27;
  int m2 = 19;
  int m3 = 25;
  unsigned long p = 2147483647;
  unsigned long m = m1 * m2 * m3;
  unsigned long r = 1;
  unsigned long bits = 150;
  unsigned long c = 2;

  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m1 << "*" << m2 << "*" << m3 << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace27 = {1426, 6176, 7126, 4751, 5701};
  const vector<long> Trace19 = {8776, 3376, 4051, 10126, 7426};
  const vector<long> Trace25 = {514, 7183, 7696, 10774, 2566, 2053};
  for (int i = 0; i < Trace27.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace27[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_16_int(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_16_int(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(256), M2(256);
  ifstream fin1, fin2;
  fin1.open("../../input/A_16_int.txt", ios::in);
  fin2.open("../../input/B_16_int.txt", ios::in);
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      fin1 >> M1[i * 16 + j];
      fin2 >> M2[i * 16 + j];
    }
  }
  cout << endl;
  cout << "The matrix A is equal to: " << endl;
  printf_matrix(M1, 16, 11, p);
  cout << endl;
  cout << "The matrix B is equal to: " << endl;
  printf_matrix(M2, 16, 11, p);
  cout << endl;

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  NTL::Vec<long> poly1, poly2;
  HELIB_NTIMER_START(Encode);
  encode_16(poly1, M1, uv);
  encode_16(poly2, M2, uv);
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  helib::Ctxt ctxt1(public_key);
  public_key.Encrypt(ctxt1, poly1);
  helib::Ctxt ctxt2(public_key);
  public_key.Encrypt(ctxt2, poly2);
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");
  helib::Ptxt<helib::BGV> plaintext_result(context);

  // Matrix addition:
  helib::Ctxt ctxt_add = ctxt1;
  HELIB_NTIMER_START(MatAdd);
  ctxt_add.addCtxt(ctxt2);
  HELIB_NTIMER_STOP(MatAdd);
  helib::printNamedTimer(std::cout, "MatAdd");

  // Decrypt:
  helib::Ptxt<helib::BGV> plaintext_result_add(context);
  HELIB_NTIMER_START(Decrypt_add);
  secret_key.Decrypt(plaintext_result_add, ctxt_add);
  HELIB_NTIMER_STOP(Decrypt_add);
  helib::printNamedTimer(std::cout, "Decrypt_add");

  // Decode:
  HELIB_NTIMER_START(Decode_add);
  vector<ZZ> result_add(256);
  decode_16_int(result_add, Base, plaintext_result_add.getPolyRepr(), uv);
  HELIB_NTIMER_STOP(Decode_add);
  helib::printNamedTimer(std::cout, "Decode_add");

  cout << endl;

  cout << "Ciphertext addition result: The matrix A + B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  printf_matrix(result_add, 16, 11, p);

  cout << endl;

  cout << "The running times (in seconds) is: " << endl;

  // Matrix Multiplication:
  helib::Ctxt ctxt = ctxt1;
  HELIB_NTIMER_START(MatMult_new);
  HEmul_16_int(ctxt, ctxt1, ctxt2, vdw, wdud);
  HELIB_NTIMER_STOP(MatMult_new);
  helib::printNamedTimer(std::cout, "MatMult_new");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt);
  secret_key.Decrypt(plaintext_result, ctxt);
  HELIB_NTIMER_STOP(Decrypt);
  helib::printNamedTimer(std::cout, "Decrypt");
  vector<ZZ> result(256);

  // Decode:
  HELIB_NTIMER_START(Decode);
  decode_16_int(result, Base, plaintext_result.getPolyRepr(), uv);
  HELIB_NTIMER_STOP(Decode);
  helib::printNamedTimer(std::cout, "Decode");

  cout << endl;

  cout << "Ciphertext addition result: The matrix A * B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  printf_matrix(result, 16, 11, p);
}

void matrix_32_int(int thread)
{
  omp_set_num_threads(thread);
  int m1 = 27;
  int m2 = 19;
  int m3 = 25;
  int M = 2;
  int d = 32;
  unsigned long p = 2147483647;
  unsigned long m = m1 * m2 * m3;
  unsigned long r = 1;
  unsigned long bits = 150;
  unsigned long c = 2;

  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m1 << "*" << m2 << "*" << m3 << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace27 = {1426, 6176, 7126, 4751, 5701};
  const vector<long> Trace19 = {8776, 3376, 4051, 10126, 7426};
  const vector<long> Trace25 = {514, 7183, 7696, 10774, 2566, 2053};
  for (int i = 0; i < Trace27.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace27[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_16_int(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_16_int(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(d * d), M2(d * d);
  ifstream fin1, fin2;
  fin1.open("../../input/A_32_int.txt", ios::in);
  fin2.open("../../input/B_32_int.txt", ios::in);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      fin1 >> M1[i * d + j];
      fin2 >> M2[i * d + j];
    }
  }
  fin1.close();
  fin2.close();
  cout << endl;
  cout << "The matrix M1 is equal to: " << endl;
  printf_matrix(M1, d, 5, p);
  cout << endl;
  cout << "The matrix M2 is equal to: " << endl;
  printf_matrix(M2, d, 5, p);
  cout << endl;

  // Block matrices
  vector<vector<vector<ZZ>>> b_M1(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  vector<vector<vector<ZZ>>> b_M2(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  for (int i = 0; i < d * d; i++)
  {
    int i0 = i / d;
    int j0 = i % d;
    b_M1[i0 / 16][j0 / 16].push_back(M1[i]);
    b_M2[i0 / 16][j0 / 16].push_back(M2[i]);
  }

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  vector<vector<NTL::Vec<long>>> poly1(M, vector<NTL::Vec<long>>(M)), poly2(M, vector<NTL::Vec<long>>(M));
  HELIB_NTIMER_START(Encode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    encode_16(poly1[i][j], b_M1[i][j], uv);
    encode_16(poly2[i][j], b_M2[i][j], uv);
  }
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  vector<vector<helib::Ctxt>> ctxt1(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  vector<vector<helib::Ctxt>> ctxt2(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    public_key.Encrypt(ctxt1[i][j], poly1[i][j]);
    public_key.Encrypt(ctxt2[i][j], poly2[i][j]);
  }
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");

  // Matrix addition:
  vector<vector<helib::Ctxt>> ctxt_add(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HELIB_NTIMER_START(MatAdd);
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < M; j++)
    {
      ctxt_add[i][j] = ctxt1[i][j];
      ctxt_add[i][j].addCtxt(ctxt2[i][j]);
    }
  }
  HELIB_NTIMER_STOP(MatAdd);
  helib::printNamedTimer(std::cout, "MatAdd");

  HELIB_NTIMER_START(Decrypt_add);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result_add(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result_add[i][j], ctxt_add[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt_add);
  helib::printNamedTimer(std::cout, "Decrypt_add");

  // Decode:
  vector<vector<vector<ZZ>>> result_add(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode_add);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_int(result_add[i][j], Base, plaintext_result_add[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode_add);
  helib::printNamedTimer(std::cout, "Decode_add");

  cout << endl;

  cout << "The matrix A + B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total_add(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total_add[i * d + j] = result_add[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  ofstream fout;
  fout.open("./M_32_add.txt", ios::out);
  printf_matrix(result_total_add, d, 11, p, fout);
  fout.close();
  cout << "   See M_32_add.txt" << endl;
  cout << endl;
  cout << "The running times (in seconds) is: " << endl;

  // Matrix Multiplication:
  HELIB_NTIMER_START(MatMult_new);
  vector<vector<helib::Ctxt>> ctxt(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HEmul_block_int(ctxt, ctxt1, ctxt2, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_new);
  helib::printNamedTimer(std::cout, "MatMult_new");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result[i][j], ctxt[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt);
  helib::printNamedTimer(std::cout, "Decrypt");

  // Decode:
  vector<vector<vector<ZZ>>> result(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_int(result[i][j], Base, plaintext_result[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode);
  helib::printNamedTimer(std::cout, "Decode");

  cout << endl;

  cout << "The matrix A * B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total[i * d + j] = result[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  fout.open("./M_32_mult.txt", ios::out);
  printf_matrix(result_total, d, 11, p, fout);
  cout << "   See M_32_mult.txt" << endl;
  fout.close();
}

void matrix_64_int(int thread)
{
  omp_set_num_threads(thread);
  int m1 = 27;
  int m2 = 19;
  int m3 = 25;
  int M = 4;
  int d = 64;
  unsigned long p = 2147483647;
  unsigned long m = m1 * m2 * m3;
  unsigned long r = 1;
  unsigned long bits = 170;
  unsigned long c = 2;

  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m1 << "*" << m2 << "*" << m3 << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace27 = {1426, 6176, 7126, 4751, 5701};
  const vector<long> Trace19 = {8776, 3376, 4051, 10126, 7426};
  const vector<long> Trace25 = {514, 7183, 7696, 10774, 2566, 2053};
  for (int i = 0; i < Trace27.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace27[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_16_int(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_16_int(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(d * d), M2(d * d);
  ifstream fin1, fin2;
  cout << "Read M1 in '../../input/M1_64.txt'" << endl;
  fin1.open("../../input/A_64_int.txt", ios::in);
  cout << "Read M2 in '../../input/M2_64.txt'" << endl;
  fin2.open("../../input/B_64_int.txt", ios::in);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      fin1 >> M1[i * d + j];
      fin2 >> M2[i * d + j];
    }
  }
  fin1.close();
  fin2.close();
  ofstream fout1, fout2;
  fout1.open("../../input/A_64_int.txt", ios::out);
  fout2.open("../../input/B_64_int.txt", ios::out);
  printf_matrix(M1, d, 5, p, fout1);
  printf_matrix(M2, d, 5, p, fout2);
  // Block matrices
  vector<vector<vector<ZZ>>> b_M1(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  vector<vector<vector<ZZ>>> b_M2(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  for (int i = 0; i < d * d; i++)
  {
    int i0 = i / d;
    int j0 = i % d;
    b_M1[i0 / 16][j0 / 16].push_back(M1[i]);
    b_M2[i0 / 16][j0 / 16].push_back(M2[i]);
  }

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  vector<vector<NTL::Vec<long>>> poly1(M, vector<NTL::Vec<long>>(M)), poly2(M, vector<NTL::Vec<long>>(M));
  HELIB_NTIMER_START(Encode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    encode_16(poly1[i][j], b_M1[i][j], uv);
    encode_16(poly2[i][j], b_M2[i][j], uv);
  }
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  vector<vector<helib::Ctxt>> ctxt1(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  vector<vector<helib::Ctxt>> ctxt2(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    public_key.Encrypt(ctxt1[i][j], poly1[i][j]);
    public_key.Encrypt(ctxt2[i][j], poly2[i][j]);
  }
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");

  // Matrix addition:
  vector<vector<helib::Ctxt>> ctxt_add(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HELIB_NTIMER_START(MatAdd);
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < M; j++)
    {
      ctxt_add[i][j] = ctxt1[i][j];
      ctxt_add[i][j].addCtxt(ctxt2[i][j]);
    }
  }
  HELIB_NTIMER_STOP(MatAdd);
  helib::printNamedTimer(std::cout, "MatAdd");

  HELIB_NTIMER_START(Decrypt_add);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result_add(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result_add[i][j], ctxt_add[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt_add);
  helib::printNamedTimer(std::cout, "Decrypt_add");

  // Decode:
  vector<vector<vector<ZZ>>> result_add(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode_add);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_int(result_add[i][j], Base, plaintext_result_add[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode_add);
  helib::printNamedTimer(std::cout, "Decode_add");

  cout << endl;

  cout << "The matrix A + B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total_add(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total_add[i * d + j] = result_add[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  ofstream fout;
  fout.open("./M_64_add.txt", ios::out);
  printf_matrix(result_total_add, d, 11, p, fout);
  fout.close();
  cout << "   See M_64_add.txt" << endl;

  cout << endl;
  cout << "The running times (in seconds) is: " << endl;

  // // Matrix Multiplication:
  HELIB_NTIMER_START(MatMult_new);
  vector<vector<helib::Ctxt>> ctxt(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  HEmul_block_int(ctxt, ctxt1, ctxt2, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_new);
  helib::printNamedTimer(std::cout, "MatMult_new");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result[i][j], ctxt[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt);
  helib::printNamedTimer(std::cout, "Decrypt");
  // Decode:
  vector<vector<vector<ZZ>>> result(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_int(result[i][j], Base, plaintext_result[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode);
  helib::printNamedTimer(std::cout, "Decode");
  cout << endl;

  // Printf result;
  cout << "The matrix A * B "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total[i * d + j] = result[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  fout.open("./M_64_mult.txt", ios::out);
  printf_matrix(result_total, d, 11, p, fout);
  cout << "   See M_64_mult.txt" << endl;
  fout.close();
}

void matrix_16_int_depth_4()
{
  int m1 = 27 * 3;
  int m2 = 17;
  int m3 = 25;
  unsigned long p = 2147483647;
  unsigned long m = m1 * m2 * m3;
  unsigned long r = 1;
  unsigned long bits = 400;
  unsigned long c = 2;

  HELIB_NTIMER_START(Context);
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m1 << "*" << m2 << "*" << m3 << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  HELIB_NTIMER_STOP(Context);
  helib::printNamedTimer(std::cout, "Context");
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace27 = {2126, 10201, 5101, 3401, 6376};
  const vector<long> Trace19 = {2026, 8101, 12151, 32401};
  const vector<long> Trace25 = {1378, 5509, 20656, 22033, 6886};
  for (int i = 0; i < Trace27.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace27[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_16_int_depth_4(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_16_int_depth_4(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(256), M2(256);
  ifstream fin1;
  fin1.open("../../input/A_16.txt", ios::in);
  for (int i = 0; i < 16; i++)
  {
    for (int j = 0; j < 16; j++)
    {
      fin1 >> M1[i * 16 + j];
    }
  }
  cout << endl;
  cout << "The matrix A is equal to: " << endl;
  printf_matrix(M1, 16, 5, p);

  // //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  NTL::Vec<long> poly1;
  HELIB_NTIMER_START(Encode);
  encode_16(poly1, M1, uv);
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  helib::Ctxt ctxt1(public_key);
  public_key.Encrypt(ctxt1, poly1);
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");
  helib::Ptxt<helib::BGV> plaintext_result(context);

  // Matrix Multiplication:
  helib::Ctxt ctxt = ctxt1;
  cout << "  Multiplication:" << endl;
  cout << "    A1 = A * A: ";
  HELIB_NTIMER_START(MatMult_A1);
  HEmul_16_int_depth_4(ctxt, ctxt1, ctxt1, vdw, wdud);
  HELIB_NTIMER_STOP(MatMult_A1);
  helib::printNamedTimer(std::cout, "MatMult_A1");
  cout << "    A2 = A1 * A1: ";
  HELIB_NTIMER_START(MatMult_A2);
  HEmul_16_int_depth_4(ctxt, ctxt, ctxt, vdw, wdud);
  HELIB_NTIMER_STOP(MatMult_A2);
  helib::printNamedTimer(std::cout, "MatMult_A2");
  cout << "    A3 = A2 * A2: ";
  HELIB_NTIMER_START(MatMult_A3);
  HEmul_16_int_depth_4(ctxt, ctxt, ctxt, vdw, wdud);
  HELIB_NTIMER_STOP(MatMult_A3);
  helib::printNamedTimer(std::cout, "MatMult_A3");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt);
  secret_key.Decrypt(plaintext_result, ctxt);
  HELIB_NTIMER_STOP(Decrypt);
  helib::printNamedTimer(std::cout, "Decrypt");
  vector<ZZ> result(256);

  // Decode:
  HELIB_NTIMER_START(Decode);
  decode_16_int_depth_4(result, Base, plaintext_result.getPolyRepr(), uv);
  HELIB_NTIMER_STOP(Decode);
  helib::printNamedTimer(std::cout, "Decode");

  cout << endl;

  cout << "The matrix A^(2^3) = A^(8) "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  printf_matrix(result, 16, 11, p);
}

void matrix_32_int_depth_4(int thread)
{
  omp_set_num_threads(thread);
  int m1 = 27 * 3;
  int m2 = 17;
  int m3 = 25;
  int M = 2;
  int d = 32;
  unsigned long p = 2147483647;
  unsigned long m = m1 * m2 * m3;
  unsigned long r = 1;
  unsigned long bits = 400;
  unsigned long c = 2;

  HELIB_NTIMER_START(Context);
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m1 << "*" << m2 << "*" << m3 << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  HELIB_NTIMER_STOP(Context);
  helib::printNamedTimer(std::cout, "Context");
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace27 = {2126, 10201, 5101, 3401, 6376};
  const vector<long> Trace19 = {2026, 8101, 12151, 32401};
  const vector<long> Trace25 = {1378, 5509, 20656, 22033, 6886};
  for (int i = 0; i < Trace27.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace27[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_16_int_depth_4(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_16_int_depth_4(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(d * d);
  ifstream fin1;
  fin1.open("../../input/A_32.txt", ios::in);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      fin1 >> M1[i * d + j];
    }
  }
  fin1.close();
  cout << endl;
  cout << "The matrix A is equal to: " << endl;
  printf_matrix(M1, d, 5, p);
  cout << endl;

  // Block matrices
  vector<vector<vector<ZZ>>> b_M1(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  for (int i = 0; i < d * d; i++)
  {
    int i0 = i / d;
    int j0 = i % d;
    b_M1[i0 / 16][j0 / 16].push_back(M1[i]);
  }

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  vector<vector<NTL::Vec<long>>> poly1(M, vector<NTL::Vec<long>>(M));
  HELIB_NTIMER_START(Encode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    encode_16(poly1[i][j], b_M1[i][j], uv);
  }
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  vector<vector<helib::Ctxt>> ctxt1(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    public_key.Encrypt(ctxt1[i][j], poly1[i][j]);
  }
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");

  // Matrix Multiplication:
  vector<vector<helib::Ctxt>> ctxt(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  cout << "  Multiplication:" << endl;
  cout << "    A1 = A * A: ";
  HELIB_NTIMER_START(MatMult_A1);
  HEmul_block_int_depth_4(ctxt, ctxt1, ctxt1, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_A1);
  helib::printNamedTimer(std::cout, "MatMult_A1");
  cout << "    A2 = A1 * A1: ";
  HELIB_NTIMER_START(MatMult_A2);
  HEmul_block_int_depth_4(ctxt, ctxt, ctxt, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_A2);
  helib::printNamedTimer(std::cout, "MatMult_A2");
  cout << "    A3 = A2 * A2: ";
  HELIB_NTIMER_START(MatMult_A3);
  HEmul_block_int_depth_4(ctxt, ctxt, ctxt, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_A3);
  helib::printNamedTimer(std::cout, "MatMult_A3");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result[i][j], ctxt[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt);
  helib::printNamedTimer(std::cout, "Decrypt");

  // Decode:
  vector<vector<vector<ZZ>>> result(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_int_depth_4(result[i][j], Base, plaintext_result[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode);
  helib::printNamedTimer(std::cout, "Decode");

  cout << endl;

  cout << "The matrix A^(2^3) = A^(8) "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total[i * d + j] = result[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  ofstream fout;
  fout.open("./M_32.txt", ios::out);
  printf_matrix(result_total, d, 11, p, fout);
  cout << "   See M_32.txt" << endl;
  fout.close();
}

void matrix_64_int_depth_4(int thread)
{
  omp_set_num_threads(thread);
  int m1 = 27 * 3;
  int m2 = 17;
  int m3 = 25;
  int M = 4;
  int d = 64;
  unsigned long p = 2147483647;
  unsigned long m = m1 * m2 * m3;
  unsigned long r = 1;
  unsigned long bits = 400;
  unsigned long c = 2;

  HELIB_NTIMER_START(Context);
  std::cout << "Initialising context object..." << std::endl;
  helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();
  std::cout << "  m = " << m1 << "*" << m2 << "*" << m3 << "\n"
            << "  number of bits of Q = " << context.bitSizeOfQ() << "\n"
            << "  plaintext modulus q = " << p << "\n"
            << "  security level = " << context.securityLevel() << std::endl;
  std::cout << "  Creating secret key..." << std::endl;
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  HELIB_NTIMER_STOP(Context);
  helib::printNamedTimer(std::cout, "Context");
  helib::SecKey secret_key(context);
  // Generate the secret key and key-switching matrices that we need
  // HeapProfilerStart("test");
  secret_key.GenSecKey();
  std::cout << "  Generating key-switching matrices..." << std::endl;
  const vector<long> Trace27 = {2126, 10201, 5101, 3401, 6376};
  const vector<long> Trace19 = {2026, 8101, 12151, 32401};
  const vector<long> Trace25 = {1378, 5509, 20656, 22033, 6886};
  for (int i = 0; i < Trace27.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace27[i], 0, 0);
  }
  for (int i = 0; i < Trace19.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace19[i], 0, 0);
  }
  for (int i = 0; i < Trace25.size(); i++)
  {
    secret_key.GenKeySWmatrix(1, Trace25[i], 0, 0);
  }
  secret_key.setKeySwitchMap();
  const helib::PubKey &public_key = secret_key;

  // load Encoded bases and decoded bases
  vector<ZZX> uv(256);
  helib::zzX vdw;
  helib::zzX wdud;
  encode_base_gen_16_int_depth_4(uv, vdw, wdud);
  vector<ZZX> Base(256);
  decode_base_gen_16_int_depth_4(Base);

  // Read matrices M1 and M2;
  vector<ZZ> M1(d * d), M2(d * d);
  ifstream fin1, fin2;
  fin1.open("../../input/A_64.txt", ios::in);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      fin1 >> M1[i * d + j];
    }
  }
  fin1.close();
  cout << endl;
  cout << "   Read A in A_64.txt." << endl;
  //printf_matrix(M1, d, 5, p);
  cout << endl;
  // Block matrices
  vector<vector<vector<ZZ>>> b_M1(M, vector<vector<ZZ>>(M, vector<ZZ>()));
  for (int i = 0; i < d * d; i++)
  {
    int i0 = i / d;
    int j0 = i % d;
    b_M1[i0 / 16][j0 / 16].push_back(M1[i]);
  }

  //-------------------------------------------------------------
  cout << "The running times (in seconds) is: " << endl;

  // Encode:
  vector<vector<NTL::Vec<long>>> poly1(M, vector<NTL::Vec<long>>(M));
  HELIB_NTIMER_START(Encode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    encode_16(poly1[i][j], b_M1[i][j], uv);
  }
  HELIB_NTIMER_STOP(Encode);
  helib::printNamedTimer(std::cout, "Encode");

  // // Encrypt:
  HELIB_NTIMER_START(Encrypt);
  vector<vector<helib::Ctxt>> ctxt1(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    public_key.Encrypt(ctxt1[i][j], poly1[i][j]);
  }
  HELIB_NTIMER_STOP(Encrypt);
  helib::printNamedTimer(std::cout, "Encrypt");

  // // Matrix Multiplication:
  HELIB_NTIMER_START(MatMult_new);
  vector<vector<helib::Ctxt>> ctxt(M, vector<helib::Ctxt>(M, helib::Ctxt(public_key)));
  cout << "  Multiplication:" << endl;
  cout << "    A1 = A * A: ";
  HELIB_NTIMER_START(MatMult_A1);
  HEmul_block_int_depth_4(ctxt, ctxt1, ctxt1, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_A1);
  helib::printNamedTimer(std::cout, "MatMult_A1");
  cout << "    A2 = A1 * A1: ";
  HELIB_NTIMER_START(MatMult_A2);
  HEmul_block_int_depth_4(ctxt, ctxt, ctxt, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_A2);
  helib::printNamedTimer(std::cout, "MatMult_A2");
  cout << "    A3 = A2 * A2: ";
  HELIB_NTIMER_START(MatMult_A3);
  HEmul_block_int_depth_4(ctxt, ctxt, ctxt, vdw, wdud, M);
  HELIB_NTIMER_STOP(MatMult_A3);
  helib::printNamedTimer(std::cout, "MatMult_A3");
  HELIB_NTIMER_STOP(MatMult_new);
  helib::printNamedTimer(std::cout, "MatMult_new");

  // Decrypt:
  HELIB_NTIMER_START(Decrypt);
  vector<vector<helib::Ptxt<helib::BGV>>> plaintext_result(M, vector<helib::Ptxt<helib::BGV>>(M, helib::Ptxt<helib::BGV>(context)));
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    secret_key.Decrypt(plaintext_result[i][j], ctxt[i][j]);
  }
  HELIB_NTIMER_STOP(Decrypt);
  helib::printNamedTimer(std::cout, "Decrypt");

  // Decode:
  vector<vector<vector<ZZ>>> result(M, vector<vector<ZZ>>(M, vector<ZZ>(256)));
  HELIB_NTIMER_START(Decode);
#pragma omp parallel for
  for (int L = 0; L < M * M; L++)
  {
    int i = L / M;
    int j = L % M;
    decode_16_int_depth_4(result[i][j], Base, plaintext_result[i][j].getPolyRepr(), uv);
  }
  HELIB_NTIMER_STOP(Decode);
  helib::printNamedTimer(std::cout, "Decode");

  cout << endl;

  cout << "The matrix A^(2^3) = A^(8) "
       << "mod " << p << " ("
       << "-" << p / 2 << "~" << p / 2 << ") "
       << " is equal to : " << endl;
  vector<ZZ> result_total(d * d);
  for (int i = 0; i < d; i++)
  {
    for (int j = 0; j < d; j++)
    {
      result_total[i * d + j] = result[i / 16][j / 16][(i % 16) * 16 + (j % 16)];
    }
  }
  ofstream fout;
  fout.open("./M_64.txt", ios::out);
  printf_matrix(result_total, d, 11, p, fout);
  cout << "   See M_64.txt" << endl;
  fout.close();
}

int main()
{
  int s;
  int t, d, thread, L;
  cout << "            ------------------------------- " << endl;
  cout << "             Momomorpic Matrix Operations" << endl;
  cout << "            ------------------------------- " << endl;
  cout << "Please chose the size of number: '0' for '-32768~32768' and '1' for '-1073741823~1073741823" << endl;
  cin >> t;
  cout << "Please input the size of matrix (16, 32 or 64):" << endl;
  cin >> d;
  if (t == 0)
  {
    if (d == 16)
      matrix_16_short_int();
    else
    {
      cout << "Please input the number of thread:" << endl;
      cin >> thread;
      if (d == 32)
        matrix_32_short_int(thread);
      else
        matrix_64_short_int(thread);
    }
  }
  else
  {
    cout << "Please input the depth (1 or 3):" << endl;
    cin >> L;
    if (L == 1)
    {
      if (d == 16)
        matrix_16_int();
      else
      {
        cout << "Please input the number of thread:" << endl;
        cin >> thread;
        if (d == 32)
          matrix_32_int(thread);
        else
          matrix_64_int(thread);
      }
    }
    else
    {
      if (d == 16)
        matrix_16_int_depth_4();
      else
      {
        cout << "Please input the number of thread:" << endl;
        cin >> thread;
        if (d == 32)
          matrix_32_int_depth_4(thread);
        else
          matrix_64_int_depth_4(thread);
      }
    }
  }
}

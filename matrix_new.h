#include <helib/helib.h>
#include <vector>
#include <fstream>

void plain_text_multiplication(std::vector<NTL::ZZ> M1, std::vector<NTL::ZZ> M2, int d);

void encode_base_gen_short_int(std::vector<NTL::ZZX> &uv, helib::zzX &vdw, helib::zzX &wdud);

void encode_16(NTL::Vec<long> &poly, std::vector<NTL::ZZ> &M, const std::vector<NTL::ZZX> &uv);

void decode_base_gen_short_int(std::vector<NTL::ZZX> &Base);

void decode_16_short_int(std::vector<NTL::ZZ> &result_vector, std::vector<NTL::ZZX> &Base, NTL::ZZX f, const std::vector<NTL::ZZX> &uv);

void HEmul_16_short_int(helib::Ctxt &result, helib::Ctxt ctxt1, helib::Ctxt ctxt2, helib::zzX &vdw, helib::zzX &wdud);

void HEmul_block_short_int(std::vector<std::vector<helib::Ctxt>> &ctxt, std::vector<std::vector<helib::Ctxt>> ctxt1, std::vector<std::vector<helib::Ctxt>> ctxt2, helib::zzX &vdw, helib::zzX &wdud, int M);

void printf_matrix(const std::vector<NTL::ZZ> &M, int n, int w, long unsigned int &p, std::ofstream &fout);

void printf_matrix(const std::vector<NTL::ZZ> &M, int n, int w, long unsigned int &p);

void encode_base_gen_16_int(std::vector<NTL::ZZX> &uv, helib::zzX &vdw, helib::zzX &wdud);

void decode_base_gen_16_int(std::vector<NTL::ZZX> &Base);

void decode_16_int(std::vector<NTL::ZZ> &result_vector, std::vector<NTL::ZZX> &Base, NTL::ZZX f, const std::vector<NTL::ZZX> &uv);

void HEmul_16_int(helib::Ctxt &result, helib::Ctxt ctxt1, helib::Ctxt ctxt2, helib::zzX &vdw, helib::zzX &wdud);

void HEmul_block_int(std::vector<std::vector<helib::Ctxt>> &result, std::vector<std::vector<helib::Ctxt>> ctxt1, std::vector<std::vector<helib::Ctxt>> ctxt2, helib::zzX &vdw, helib::zzX &wdud, int M);

void encode_base_gen_16_int_depth_4(std::vector<NTL::ZZX> &uv, helib::zzX &vdw, helib::zzX &wdud);

void decode_base_gen_16_int_depth_4(std::vector<NTL::ZZX> &Base);

void decode_16_int_depth_4(std::vector<NTL::ZZ> &result_vector, std::vector<NTL::ZZX> &Base, NTL::ZZX f, const std::vector<NTL::ZZX> &uv);

void HEmul_16_int_depth_4(helib::Ctxt &result, helib::Ctxt ctxt1, helib::Ctxt ctxt2, helib::zzX &vdw, helib::zzX &wdud);

void HEmul_block_int_depth_4(std::vector<std::vector<helib::Ctxt>> &ctxt, std::vector<std::vector<helib::Ctxt>> ctxt1, std::vector<std::vector<helib::Ctxt>> ctxt2, helib::zzX &vdw, helib::zzX &wdud, int M);
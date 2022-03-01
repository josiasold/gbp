#ifndef LA_TOOLS_HPP_
#define LA_TOOLS_HPP_

#include <algorithm> //min
#include <iostream>

#include <NTL/mat_GF2.h> 

#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

inline void gf2_syndrome(xt::xarray<int> &s,const xt::xarray<int> &y,const xt::xarray<int> &H)
{
  s = xt::linalg::dot(H, y) % 2;
};

inline xt::xarray<int> gf2_syndrome(const xt::xarray<int> &y,const xt::xarray<int> &H)
{
  xt::xarray<int> s = xt::linalg::dot(H, y) % 2;
  return s;
};

inline void gf2_rank(const xt::xarray<int> &Matrix, int n_c, int n_q, int &r)
{
  NTL::Mat<NTL::GF2> M;
  M.SetDims(n_c,n_q);
  for (long i = 0; i < n_c; i++) 
  {
    for (long j = 0; j < n_q; j++) 
    {
      M[i][j] = Matrix(i,j);
    }
  }
  r = NTL::gauss(M);
};
inline int gf2_rank(const xt::xarray<int> &Matrix, int n_c, int n_q)
{
  NTL::Mat<NTL::GF2> M;
  M.SetDims(n_c,n_q);
  for (long i = 0; i < n_c; i++) 
  {
    for (long j = 0; j < n_q; j++) 
    {
      M[i][j] = Matrix(i,j);
    }
  }
  int r = NTL::gauss(M);
  return r;
};

inline bool gf2_isEquiv(const xt::xarray<int> &e, const xt::xarray<int> &H, int n_c, int n_q)
{
  int rank_init = gf2_rank(H, n_c, n_q);

  xt::xarray<int> H_larger = xt::zeros_like(H);
  H_larger.resize({static_cast<unsigned long>(n_c + 1), static_cast<unsigned long>(n_q)});
  xt::view(H_larger, xt::range(0, n_c), xt::all()) = H;
  xt::row(H_larger, n_c) = e;

  int rank_after = gf2_rank(H_larger, n_c + 1, n_q);
  if (rank_after == rank_init)
    return true;
  else
    return false;
};

inline void gf4_syndrome(xt::xarray<int> &s, const xt::xarray<int> &y, const xt::xarray<int> &H){

  static const int mul_table[16] = {0, 0, 0, 0, 0, 1, 2, 3, 0, 2, 3, 1, 0, 3, 1, 2};
  static const int conj_table[4] = {0, 1, 3, 2};
  static const int trace_table[4] = {0, 0, 1, 1};
  static const int stride = 4;

  int n_c = H.shape(0);
  int n_q = H.shape(1);

  for (size_t c = 0; c < n_c; c++) 
  {
    int s_value = 0;
    for (size_t q = 0; q < n_q; q++) 
    {
        s_value ^= mul_table[H(c,q) * stride + conj_table[y(q)]];
    }
    s(c) = trace_table[s_value];
  }
};

inline xt::xarray<int> gf4_syndrome(const xt::xarray<int> &y, const xt::xarray<int> &H)
{
  static const int mul_table[16] = {0, 0, 0, 0, 0, 1, 2, 3, 0, 2, 3, 1, 0, 3, 1, 2};
  static const int conj_table[4] = {0, 1, 3, 2};
  static const int trace_table[4] = {0, 0, 1, 1};
  static const int stride = 4;

  int n_c = H.shape(0);
  int n_q = H.shape(1);

  xt::xarray<int> s = xt::zeros<int>({n_c});

  for (size_t c = 0; c < n_c; c++) 
  {
    int s_value = 0;
    for (size_t q = 0; q < n_q; q++) 
    {
        s_value ^= mul_table[H(c,q) * stride + conj_table[y(q)]];
    }
    s(c) = trace_table[s_value];
  }
  return s;
};

inline void gf4_rank(const xt::xarray<int> &Matrix, int n_c, int n_q, int &r)
{
  NTL::Mat<NTL::GF2> M_X;
  M_X.SetDims(n_c,n_q);

  NTL::Mat<NTL::GF2> M_Z;
  M_Z.SetDims(n_c,n_q);

  for (long i = 0; i < n_c; i++) 
  {
    for (long j = 0; j < n_q; j++) 
    {
      if (Matrix(i,j) == 1) 
      {
        M_X[i][j] = 1;
        M_Z[i][j] = 0;
      }
      else if (Matrix(i,j) == 2) 
      {
        M_X[i][j] = 0;
        M_Z[i][j] = 1;
      }
      else if (Matrix(i,j) == 3) 
      {
        M_X[i][j] = 1;
        M_Z[i][j] = 1;
      }
      else
      {
        M_X[i][j] = 0;
        M_Z[i][j] = 0;
      }
    }
  }

  int r_X = NTL::gauss(M_X);
  int r_Z = NTL::gauss(M_Z);

  r = r_X + r_Z;
};

inline int gf4_rank(const xt::xarray<int> &Matrix, int n_c, int n_q)
{
  NTL::Mat<NTL::GF2> M_X;
  M_X.SetDims(n_c,n_q);

  NTL::Mat<NTL::GF2> M_Z;
  M_Z.SetDims(n_c,n_q);

  for (long i = 0; i < n_c; i++) 
  {
    for (long j = 0; j < n_q; j++) 
    {
      if (Matrix(i,j) == 1) 
      {
        M_X[i][j] = 1;
        M_Z[i][j] = 0;
      }
      else if (Matrix(i,j) == 2) 
      {
        M_X[i][j] = 0;
        M_Z[i][j] = 1;
      }
      else if (Matrix(i,j) == 3) 
      {
        M_X[i][j] = 1;
        M_Z[i][j] = 1;
      }
      else
      {
        M_X[i][j] = 0;
        M_Z[i][j] = 0;
      }
    }
  }

  int r_X = NTL::gauss(M_X);
  int r_Z = NTL::gauss(M_Z);

  int r = r_X + r_Z;
  return r;
};

inline bool gf4_isEquiv(const xt::xarray<int> &e, const xt::xarray<int> &H, int n_c, int n_q)
{
  int rank_init = gf4_rank(H, n_c, n_q);

  xt::xarray<int> H_larger = xt::zeros_like(H);
  H_larger.resize({static_cast<unsigned long>(n_c + 1), static_cast<unsigned long>(n_q)});
  xt::view(H_larger, xt::range(0, n_c), xt::all()) = H;
  xt::row(H_larger, n_c) = e;

  int rank_after = gf4_rank(H_larger, n_c + 1, n_q);
  if (rank_after == rank_init)
    return true;
  else
    return false;
};

inline int gf4_mul(int a, int b)
{
  static const int mul_table[16] = {0, 0, 0, 0, 0, 1, 2, 3, 0, 2, 3, 1, 0, 3, 1, 2};
  return mul_table[a*4+b];
};

inline int gf4_conj(int a)
{
  static const int conj_table[4] = {0, 1, 3, 2};
  return conj_table[a];
};

inline int hamming_weight(const xt::xarray<int> &x)
{
  int hw = 0;
  for (auto it = x.begin(); it != x.end(); ++it)
  {
    if (*it != 0) hw++;
  }
  return hw;
};

inline xt::xarray<int> get_x(const xt::xarray<int> &y)
{
  xt::xarray<int> x = xt::zeros_like(y);
  for (int i = 0; i<y.size(); i++)
  {
    if ((y(i) == 1) || (y(i) == 3))
    {
      x(i) = 1;
    }
  }
  return x;
};

inline xt::xarray<int> get_z(const xt::xarray<int> &y)
{
  xt::xarray<int> z = xt::zeros_like(y);
  for (int i = 0; i<y.size(); i++)
  {
    if ((y(i) == 2) || (y(i) == 3))
    {
      z(i) = 1;
    }
  }
  return z;
};

#endif
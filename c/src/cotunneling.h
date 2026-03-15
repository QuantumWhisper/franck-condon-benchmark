#ifndef FC_COTUNNELING_H
#define FC_COTUNNELING_H

#include "fc_matrix.h"

void m_sumMMr(int N,
              int q1,
              const int *q2_vec,
              int nq2,
              double lambda,
              double muL,
              double muR,
              double vmode,
              double epsilond,
              double T,
              FCCache *fc,
              double *out);

void m_sumMMr11(int N,
                int q1,
                const int *q2_vec,
                int nq2,
                double lambda,
                double muL,
                double muR,
                double vmode,
                double epsilond,
                double T,
                FCCache *fc,
                double *out);

double m_sumMMMMrs(int N,
                   int q1,
                   int q2,
                   double lambda,
                   double muL,
                   double muR,
                   double vmode,
                   double epsilond,
                   double T,
                   FCCache *fc);

double m_sumMMMMrs11(int N,
                     int q1,
                     int q2,
                     double lambda,
                     double muL,
                     double muR,
                     double vmode,
                     double epsilond,
                     double T,
                     FCCache *fc);

void sumMMr(int q1,
            const int *q2_vec,
            int nq2,
            double lambda,
            double muL,
            double muR,
            double vmode,
            double epsilond,
            double T,
            FCCache *fc,
            double *out);

void sumMMr11(int q1,
              const int *q2_vec,
              int nq2,
              double lambda,
              double muL,
              double muR,
              double vmode,
              double epsilond,
              double T,
              FCCache *fc,
              double *out);

void sumMMMMrs(int q1,
               const int *q2_vec,
               int nq2,
               double lambda,
               double muL,
               double muR,
               double vmode,
               double epsilond,
               double T,
               FCCache *fc,
               double *out);

void sumMMMMrs11(int q1,
                 const int *q2_vec,
                 int nq2,
                 double lambda,
                 double muL,
                 double muR,
                 double vmode,
                 double epsilond,
                 double T,
                 FCCache *fc,
                 double *out);

#endif

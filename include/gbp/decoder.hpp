#ifndef DECODER_HPP_
#define DECODER_HPP_

#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp> // for xt::row
#include <xtensor/xnpy.hpp> // for xt::dump_npy

#include "region_graph.hpp"
#include "la_tools.hpp"
#include "io_tools.hpp"
#include "timing.hpp"

namespace gbp
{
        class Decoder
        {
        private:
                /* data */
                const int p_max_iterations;
                int p_took_iterations;
                const int p_max_repetitions;
                int p_took_repetitions;
                long double p_w;

                xt::xarray<int> p_H;
                const xt::xarray<int> p_H_X;
                const xt::xarray<int> p_H_Z;

                int p_n_qubits;
                int p_n_checks;
                int p_n_checks_X;
                int p_n_checks_Z;

                gbp::Regiongraph p_regionGraph_X;
                gbp::Regiongraph p_regionGraph_Z;

                bool p_save_history;
                xt::xarray<int> p_error_guess_history;
                xt::xarray<int> p_syndrome_history;
                xt::xarray<int> p_incompatibility_score_history;
                xt::xarray<long double> p_tdqs_history;

                int p_print_detail;

        public:
                Decoder(xt::xarray<int> H_X, xt::xarray<int> H_Z,const int max_iterations,const int max_repetitions,const long double w, bool save_history, int print_detail): p_regionGraph_X(H_Z), p_regionGraph_Z(H_X), p_max_iterations(max_iterations), p_max_repetitions(max_repetitions), p_w(w), p_H_X(H_X), p_H_Z(H_Z), p_save_history(save_history), p_print_detail(print_detail)
                {
                        /// p_regionGraph_X is used for decoding X-errors --> has to use H_Z = Z-stabilizers
                        p_n_qubits = H_X.shape(1);
                        p_n_checks_X = H_X.shape(0);
                        p_n_checks_Z = H_Z.shape(0);
                        p_n_checks = p_n_checks_X + p_n_checks_Z;
                        p_H = xt::zeros<int>({p_n_checks, p_n_qubits});
                        auto p_H_X_view = xt::view(p_H, xt::range(0, p_n_checks_X), xt::all());
                        auto p_H_Z_view = xt::view(p_H, xt::range(p_n_checks_X, p_n_checks), xt::all());

                        p_H_X_view = H_X;
                        p_H_Z_view = 2 * H_Z;

                        if (p_save_history)
                        {
                                int length = p_max_repetitions* (p_max_iterations+1);
                                p_error_guess_history = -1*xt::ones<int>({length, p_n_qubits});
                                p_syndrome_history = -1*xt::ones<int>({length, p_n_checks});
                                p_tdqs_history = -1*xt::ones<int>({length,3});
                                p_incompatibility_score_history = xt::zeros<int>({length});
                        }

                        p_regionGraph_X.set_print_detail(print_detail);
                        p_regionGraph_Z.set_print_detail(print_detail);
                } ;

                xt::xarray<int> decode(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome_0);
                void dump_history(std::string PATH)
                {
                        std::string SYNDROME_PATH = PATH + "/syndromes.npy";
                        std::string ERRORGUESS_PATH = PATH + "/errorguesses.npy";
                        std::string TDQS_PATH = PATH + "/tdqs.npy";
                        std::string INCOMPATIBILITYSCORE_PATH = PATH + "/incompatibility.npy";
                        xt::dump_npy(SYNDROME_PATH,p_syndrome_history);
                        xt::dump_npy(ERRORGUESS_PATH,p_error_guess_history);
                        xt::dump_npy(TDQS_PATH,p_tdqs_history);
                        xt::dump_npy(INCOMPATIBILITYSCORE_PATH,p_incompatibility_score_history);
                }

                int took_repetitions(){return p_took_repetitions;};
                int took_iterations(){return p_took_iterations;};
        };
        
} // end of namespace gbp

#endif
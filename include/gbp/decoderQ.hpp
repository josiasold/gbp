#ifndef DECODERQ_HPP_
#define DECODERQ_HPP_

#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp> // for xt::row
#include <xtensor/xnpy.hpp> // for xt::dump_npy

#include "region_graphQ.hpp"
#include "tanner_graph.hpp"
#include "la_tools.hpp"
#include "io_tools.hpp"
#include "timing.hpp"
#include "properties.hpp"
#include "error_channels.hpp"

namespace gbp
{
        class DecoderQ
        {
        private:
                /* data */
                int p_took_iterations;
                int p_took_repetitions;

                xt::xarray<int> p_H;
                const xt::xarray<int> p_H_X;
                const xt::xarray<int> p_H_Z;

                int p_n_qubits;
                int p_n_checks;
                int p_n_checks_X;
                int p_n_checks_Z;

                gbp::QRegiongraph p_regionGraph;

                xt::xarray<int> p_error_guess_history;
                xt::xarray<int> p_syndrome_history;
                xt::xarray<int> p_incompatibility_score_history;
                xt::xarray<long double> p_tdqs_history;


                struct Properties 
                {
                        int max_iterations;
                        int max_repetitions;
                        size_t repeat_p0_strategy;
                        size_t repeat_stopping_criterion;
                        int hard_decision_method;
                        long double damping;
                        bool normalize;
                        bool save_history;
                        int verbose;
                } properties;


        public:
                DecoderQ(xt::xarray<int> H, xt::xarray<int> H_X, xt::xarray<int> H_Z): p_regionGraph(H), p_H(H), p_H_X(H_X), p_H_Z(H_Z)
                {
                        p_n_qubits = H_X.shape(1);
                        p_n_checks_X = H_X.shape(0);
                        p_n_checks_Z = H_Z.shape(0);
                        p_n_checks = p_n_checks_X + p_n_checks_Z;

                        setProperties({});
                } ;

                DecoderQ(xt::xarray<int> H, xt::xarray<int> H_X, xt::xarray<int> H_Z, const gbp::PropertySet &opts): p_regionGraph(H), p_H(H), p_H_X(H_X), p_H_Z(H_Z)
                {
                        p_n_qubits = H_X.shape(1);
                        p_n_checks_X = H_X.shape(0);
                        p_n_checks_Z = H_Z.shape(0);
                        p_n_checks = p_n_checks_X + p_n_checks_Z;

                        setProperties(opts);
                } ;

                xt::xarray<int> decode(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome_0, std::string channel="depolarizing");

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

                void setProperties( const gbp::PropertySet &opts );
                gbp::PropertySet getProperties() const;
                std::string printProperties() const;
        };

} // end of namespace gbp

#endif
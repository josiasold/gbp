#ifndef TANNER_GRAPH_HPP_
#define TANNER_GRAPH_HPP_

#include <lemon/list_graph.h>
#include <lemon/adaptors.h> // for mapFill
#include <xtensor/xarray.hpp> //for xt::xarray, xt::sum
#include <xtensor/xadapt.hpp> //for xt::xadapt
#include <xtensor/xsort.hpp> //for xt::xunique
#include <xtensor/xview.hpp> // for xt::row
#include <xtensor/xio.hpp> // for <<

#include "varset.hpp"
#include "factor.hpp"
#include "la_tools.hpp"
#include "io_tools.hpp"

namespace gbp
{

        class Tannergraph
        {
        private:
                /* data */
                xt::xarray<int> p_H;
                int p_n_checks;
                int p_n_qubits;

                int p_print_detail;

                int p_incompatibility_score;

                lemon::ListGraph p_tannergraph;

                lemon::ListGraph::NodeMap< int > p_node_type;
                lemon::ListGraph::NodeMap< int > p_node_index;
                lemon::ListGraph::NodeMap< int > p_counting_number;
                lemon::ListGraph::NodeMap< gbp::VarSet > p_node_qubits;
                lemon::ListGraph::NodeMap< gbp::TFactor<long double> > p_node_belief;
                lemon::ListGraph::NodeMap< gbp::TFactor<long double> > p_node_initial_belief;

                xt::xarray<int> p_hard_decision;
                xt::xarray<int> p_syndrome;
                xt::xarray<long double> p_initial_error_probabilities;

                lemon::ListGraph::NodeMap< gbp::TFactor<long double> > p_node_localfactors_lnsum; // for tdqs
                lemon::ListGraph::NodeMap< xt::xarray<long double> > p_node_tdqs; // thermodynamic quantities
                xt::xarray<long double> p_tdqs; // overall thermodynamic quantities

                lemon::ListGraph::EdgeMap< int > p_edge_type;
                lemon::ListGraph::EdgeMap< int > p_edge_id;

                lemon::ListGraph::EdgeMap< xt::xarray<int> > p_edge_H_indices;

                lemon::ListGraph::EdgeMap< long double > p_m_cq;
                lemon::ListGraph::EdgeMap< long double > p_m_qc;

                lemon::ListGraph::EdgeMap< xt::xarray<long double> > p_r;

                lemon::ListGraph::EdgeMap< bool > p_edge_converged;

                void construct_tannergraph();
                void filter_parity(int parity, int check_type, gbp::TFactor<long double> &belief);
                const xt::xarray<int> p_parity_map = {{0,0,0,0},{0,0,1,1},{0,1,0,1},{0,1,1,0}};
        public:
                Tannergraph(xt::xarray<int> t_H) : p_H(t_H), p_node_type(p_tannergraph), p_node_index(p_tannergraph), p_counting_number(p_tannergraph), p_node_qubits(p_tannergraph), p_node_belief(p_tannergraph), p_node_initial_belief(p_tannergraph), p_edge_type(p_tannergraph), p_edge_id(p_tannergraph), p_m_cq(p_tannergraph), p_m_qc(p_tannergraph), p_r(p_tannergraph), p_node_localfactors_lnsum(p_tannergraph), p_node_tdqs(p_tannergraph), p_edge_H_indices(p_tannergraph), p_edge_converged(p_tannergraph)
                {
                        p_n_checks = t_H.shape(0);
                        p_n_qubits = t_H.shape(1);
                        construct_tannergraph();

                        p_hard_decision = xt::zeros<int>({p_n_qubits});
                        p_syndrome = xt::zeros<int>({p_n_checks});
                        lemon::mapFill(p_tannergraph,p_r,xt::zeros<long double>({2}));
                        p_initial_error_probabilities = xt::zeros<long double>({4});

                        p_tdqs = xt::zeros<long double>({3});
                        lemon::mapFill(p_tannergraph,p_node_tdqs,p_tdqs);
                        lemon::mapFill(p_tannergraph,p_edge_converged,0);
                        // lemon::mapFill(p_tannergraph,p_edge_H_indices,xt::zeros<int>({2}));
                };

                int n_qubits() {return p_n_qubits;};
                int n_checks() {return p_n_checks;};

                void initialize_data(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome);

                void update_messages(const long double w);
                void update_beliefs(int harddecision_method = 1);
                
                xt::xarray<int> get_hard_decision();
                xt::xarray<long double> get_tdqs(bool update = true); // update thermodynamic quantities
                
                void set_print_detail(int print_detail){p_print_detail = print_detail;};

                void get_subH(xt::xarray<int>& H_sub, xt::xarray<int> &checks_sub, xt::xarray<int> &qubits_sub, int pauli_type);
        };

} // end of namespace gbp

#endif
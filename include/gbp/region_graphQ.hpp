#ifndef REGION_GRAPH_HPP_
#define REGION_GRAPH_HPP_

#include <set> //for std::set, set.contains
#include <vector>

#include <lemon/list_graph.h>
#include <lemon/planarity.h> // for PlanarDrawing
#include <lemon/adaptors.h> // for Undirector
#include <lemon/graph_to_eps.h>

#include <xtensor/xarray.hpp> //for xt::xarray, xt::sum
#include <xtensor/xview.hpp> // for xt::row
#include <xtensor/xio.hpp> // for <<

#include "varset.hpp"
#include "factor.hpp"
#include "la_tools.hpp"
#include "io_tools.hpp"

typedef gbp::TFactor<long double> Factor;

namespace gbp
{
        class QRegiongraph
        {
        private:
                /* data */
                xt::xarray<int> p_H;
                xt::xarray<int> p_check_types;
                int p_n_checks;
                int p_n_qubits;

                int p_print_detail;

                int p_incompatibility_score;

                lemon::ListDigraph p_region_graph;

                lemon::ListDigraph::NodeMap< int > p_region_layer;
                lemon::ListDigraph::NodeMap< long double > p_region_counting_number;
                lemon::ListDigraph::NodeMap< int > p_region_check;
                lemon::ListDigraph::NodeMap< gbp::VarSet > p_region_qubits;
                lemon::ListDigraph::NodeMap< Factor > p_region_factor;
                
                lemon::ListDigraph::NodeMap< Factor > p_region_belief;
                xt::xarray<int> p_hard_decision;

                lemon::ListDigraph::NodeMap< xt::xarray<long double> > p_region_tdqs; // thermodynamic quantities
                xt::xarray<long double> p_tdqs; // overall thermodynamic quantities

                lemon::ListDigraph::NodeMap< Factor > p_region_localfactors_product;
                lemon::ListDigraph::NodeMap< Factor > p_region_localfactors_lnsum;

                lemon::ListDigraph::ArcMap< Factor > p_message;
                lemon::ListDigraph::ArcMap< bool > p_message_converged;

                lemon::ListDigraph::NodeMap< std::vector< lemon::ListDigraph::Arc > > p_belief_constituency_list;
                lemon::ListDigraph::ArcMap< std::vector< lemon::ListDigraph::Arc > > p_message_constituency_list;

                void construct_regiongraph();

                lemon::ListDigraph::NodeMap< std::vector< Factor > > p_region_parity_filters;
                void filter_parity(int parity, int check_type, Factor &belief);
                void filter_parity(int parity, lemon::ListDigraph::Node &region);
        public:
                QRegiongraph(xt::xarray<int> t_H) : p_H(t_H),p_region_layer(p_region_graph), p_region_counting_number(p_region_graph), p_region_check(p_region_graph), p_region_qubits(p_region_graph), p_region_factor(p_region_graph), p_region_belief(p_region_graph), p_region_tdqs(p_region_graph), p_region_localfactors_product(p_region_graph), p_region_localfactors_lnsum(p_region_graph),p_message(p_region_graph), p_message_converged(p_region_graph), p_belief_constituency_list(p_region_graph), p_message_constituency_list(p_region_graph), p_region_parity_filters(p_region_graph)
                {
                        p_n_checks = t_H.shape(0);
                        p_n_qubits = t_H.shape(1);
                        p_check_types = xt::ones<int>({p_n_checks});
                        construct_regiongraph();

                        p_tdqs = xt::zeros<long double>({3});
                        lemon::mapFill(p_region_graph,p_region_tdqs,p_tdqs);
                        p_hard_decision = xt::zeros<int>({p_n_qubits});
                };

                int n_qubits() {return p_n_qubits;};
                int n_checks() {return p_n_checks;};

                void save_regiongraph(std::string suffix);
                void print_regiongraph();
                void initialize_data(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome);

                void update_messages(const long double damping, bool normalize = true);
                void update_beliefs(bool normalize = true);
                xt::xarray<long double> get_tdqs(bool update = true, bool normalize = true); // update thermodynamic quantities
                xt::xarray<int> get_hard_decision(int method = 1);
                int get_incompatibility_score(){return p_incompatibility_score;};
                
                void set_print_detail(int print_detail){p_print_detail = print_detail;};
        };


} // end of namespace gbp

#endif
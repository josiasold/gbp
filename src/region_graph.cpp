#include <gbp/region_graph.hpp>

void gbp::Regiongraph::construct_regiongraph()
{
        // superregions: for every check one region with a factor dependent on all check vars
        for (int c = 0; c < p_n_checks; c++)
        {
                std::vector<gbp::Var> check_vars;
                size_t check_degree = xt::sum(xt::row(p_H, c))();
                check_vars.reserve(check_degree);
                for (size_t q = 0; q < p_n_qubits; q++)
                {
                        if (p_H(c, q) != 0)
                        {
                                check_vars.push_back(gbp::Var(q, 2));
                        }
                }
                lemon::ListDigraph::Node new_region = p_region_graph.addNode();
                p_region_checks[new_region] = {c};
                p_region_qubits[new_region] = gbp::VarSet(check_vars.begin(), check_vars.end(), check_vars.size());
                p_region_factor[new_region] = gbp::TFactor<long double>(p_region_qubits[new_region]);
                p_region_layer[new_region] = 0;
        }

        // small regions: take intersections of all regions
        /// TODO: make more efficient by exploiting structure of H
        gbp::VarSet used_vars;
        for (lemon::ListDigraph::NodeIt superregion(p_region_graph); superregion != lemon::INVALID; ++superregion)
        {
                for (lemon::ListDigraph::NodeIt other_superregion(p_region_graph); other_superregion != lemon::INVALID; ++other_superregion)
                {
                        if (superregion != other_superregion)
                        {
                                gbp::VarSet intersection = p_region_qubits[superregion] & p_region_qubits[other_superregion];
                                if (!intersection.empty())
                                {
                                        if (!used_vars.contains(intersection.front()))
                                        {
                                                used_vars.insert(intersection.front());
                                                lemon::ListDigraph::Node new_region = p_region_graph.addNode();
                                                p_region_checks[new_region] = {-1};
                                                p_region_qubits[new_region] = intersection;
                                                p_region_factor[new_region] = gbp::TFactor<long double>(p_region_qubits[new_region]);
                                                p_region_layer[new_region] = 1;
                                        }
                                }
                        }
                }
        }

        // draw edges
        /// TODO: make more efficient by drawing edges directly when construction the new nodes
        for (lemon::ListDigraph::NodeIt r1(p_region_graph); r1 != lemon::INVALID; ++r1)
        {
                if (p_region_layer[r1] == 0)
                {
                        for (lemon::ListDigraph::NodeIt r2(p_region_graph); r2 != lemon::INVALID; ++r2)
                        {
                                if (p_region_layer[r2] == 1)
                                {
                                        gbp::VarSet intersection = p_region_qubits[r1] & p_region_qubits[r2];
                                        if (!intersection.empty())
                                        {
                                                lemon::ListDigraph::Arc new_edge = p_region_graph.addArc(r1, r2);
                                        }
                                }
                        }
                }
        }

        // fill belief_constituency_list
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                // messages from parents
                for (lemon::ListDigraph::InArcIt edge_from_parent(p_region_graph,region); edge_from_parent != lemon::INVALID; ++edge_from_parent)
                {
                        p_belief_constituency_list[region].push_back(edge_from_parent);
                }
                 // messages into children (surface, xzzx: no other descendants) from other parents
                for (lemon::ListDigraph::OutArcIt edge_to_child(p_region_graph,region); edge_to_child != lemon::INVALID; ++edge_to_child)
                {
                        for (lemon::ListDigraph::InArcIt edge_from_parent_of_child(p_region_graph,p_region_graph.target(edge_to_child)); edge_from_parent_of_child != lemon::INVALID; ++edge_from_parent_of_child)
                        {
                                if (p_region_graph.id(region) != p_region_graph.id(p_region_graph.source(edge_from_parent_of_child)))
                                {
                                        p_belief_constituency_list[region].push_back(edge_from_parent_of_child);
                                }
                        }
                }

        }

        // fill message_constituency_list
        for (lemon::ListDigraph::ArcIt edge(p_region_graph); edge != lemon::INVALID; ++edge)
        {
                std::set<int> ids = {};
                ids.insert(p_region_graph.id(edge));
                for (auto e : p_belief_constituency_list[p_region_graph.target(edge)])
                {
                        auto search = ids.find(p_region_graph.id(e));
                        if (search == ids.end())
                        {
                                ids.insert(p_region_graph.id(e));
                        }
                }
                for (auto e : p_belief_constituency_list[p_region_graph.source(edge)])
                {
                        auto search = ids.find(p_region_graph.id(e));
                        if (search == ids.end())
                        {
                                p_message_constituency_list[edge].emplace_back(e);
                                ids.insert(p_region_graph.id(e));
                        }
                }
        }

        // calculate counting numbers
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                if (p_region_layer[region] == 0)
                {
                        p_region_counting_number[region] = 1.0;
                }
                else
                {
                        p_region_counting_number[region] = -1.0;
                }
        }
}

void gbp::Regiongraph::construct_regiongraph_cvm(size_t n_checks_per_superregion, size_t cvm_type)
{
        xt::random::seed(time(NULL));
        // superregions: for every check one region with a factor dependent on all check vars
        xt::xarray<int> check_list = xt::arange<int>(p_n_checks);
        if (cvm_type == 0)
        {
                xt::random::shuffle(check_list);
        }

        std::vector<int> used_checks;
        int n_r0 =  ceil((double)p_n_checks / (double) n_checks_per_superregion);
        int c_p = 0;
        for (int i_r0 = 0; i_r0 < n_r0; i_r0++) // loop through superregions
        {
                lemon::ListDigraph::Node new_region = p_region_graph.addNode();
                // std::cout << "region " << p_region_graph.id(new_region) << std::endl;
                
                std::vector<gbp::Var> check_vars;

                std::vector<int> checks_in_r0;
                std::set<int> qubits_in_r0;

                for (int c = c_p; c < c_p + n_checks_per_superregion; c++) // add checks of check list
                {
                        int c0 = check_list(c);
                        checks_in_r0.push_back(c0);
                        // std::vector<gbp::Var> check_vars;
                        size_t check_degree = xt::sum(xt::row(p_H, c0))();
                        // check_vars.reserve(check_degree);
                        for (size_t q = 0; q < p_n_qubits; q++)
                        {
                                if (p_H(c0, q) != 0)
                                {
                                        check_vars.push_back(gbp::Var(q, 2));
                                }
                        }

                // xt::xarray<int> row = xt::row(H,c0);
                // xt::xarray<int> qubits = xt::flatten(xt::from_indices(xt::argwhere(row > 0)));
                // for (auto q = qubits.begin(); q != qubits.end(); ++q)
                // {
                        // qubits_in_r0.insert(*q);
                // }
                        // std::cout << "|   check " << c0 << " check_vars: " << container_to_string(check_vars) << std::endl;
                        if (c == p_n_checks-1) break;
                }
                c_p+=n_checks_per_superregion;
                auto last = std::unique(check_vars.begin(), check_vars.end());
                check_vars.erase(last, check_vars.end());
                // std::cout << "|   check_vars: " << check_vars << std::endl;
                p_region_qubits[new_region] = gbp::VarSet(check_vars.begin(), check_vars.end(), check_vars.size());
                p_region_checks[new_region] = xt::adapt(checks_in_r0);
                p_region_factor[new_region] = gbp::TFactor<long double>(p_region_qubits[new_region]);
                p_region_layer[new_region] = 0;
                // std::cout << " : \n\t p_region_qubits[new_region] = " << p_region_qubits[new_region] << std::endl;
                // std::cout << " : \n\t p_region_factor[new_region] = " << p_region_factor[new_region] << std::endl;
        }

        // for (size_t c = 0; c < p_n_checks; c++)
        // {
        //         std::vector<gbp::Var> check_vars;
        //         size_t check_degree = xt::sum(xt::row(p_H, c))();
        //         check_vars.reserve(check_degree);
        //         for (size_t q = 0; q < p_n_qubits; q++)
        //         {
        //                 if (p_H(c, q) != 0)
        //                 {
        //                         check_vars.push_back(gbp::Var(q, 2));
        //                 }
        //         }
        //         lemon::ListDigraph::Node new_region = p_region_graph.addNode();
        //         p_region_check[new_region] = c;
        //         p_region_qubits[new_region] = gbp::VarSet(check_vars.begin(), check_vars.end(), check_vars.size());
        //         p_region_factor[new_region] = gbp::TFactor<long double>(p_region_qubits[new_region]);
        //         p_region_layer[new_region] = 0;
        // }

        // small regions: take intersections of all regions
        size_t current_layer = 1;
        bool done = false;
        gbp::VarSet used_vars;

        while (done == false)
        {
                gbp::SmallSet<gbp::VarSet> intersections = find_intersections(current_layer-1);
                std::cout << "intersections = " << intersections << std::endl;
                if (intersections.size() > 0)
                {
                        for (auto it = intersections.begin(); it != intersections.end(); ++it)
                        {
                                lemon::ListDigraph::Node new_region = p_region_graph.addNode();
                                p_region_checks[new_region] = {-1};
                                p_region_qubits[new_region] = *it;
                                p_region_layer[new_region] = current_layer;
                                p_region_factor[new_region] = gbp::TFactor<long double>(p_region_qubits[new_region]);
                                // if ((*it).size() == 1)
                                // {
                                // single_qubits_used(it->at(0)) = -1;
                                // }
                        }
                        current_layer++;
                }
                else
                {
                        done = true;
                }
        }


        // for (lemon::ListDigraph::NodeIt superregion(p_region_graph); superregion != lemon::INVALID; ++superregion)
        // {
        //         for (lemon::ListDigraph::NodeIt other_superregion(p_region_graph); other_superregion != lemon::INVALID; ++other_superregion)
        //         {
        //                 if (superregion != other_superregion)
        //                 {
        //                         gbp::VarSet intersection = p_region_qubits[superregion] & p_region_qubits[other_superregion];
        //                         if (!intersection.empty())
        //                         {
        //                                 if (!used_vars.contains(intersection.front()))
        //                                 {
        //                                         used_vars.insert(intersection.front());
        //                                         lemon::ListDigraph::Node new_region = p_region_graph.addNode();
        //                                         p_region_check[new_region] = -1;
        //                                         p_region_qubits[new_region] = intersection;
        //                                         p_region_factor[new_region] = gbp::TFactor<long double>(p_region_qubits[new_region]);
        //                                         p_region_layer[new_region] = 1;
        //                                 }
        //                         }
        //                 }
        //         }
        // }

        // draw edges
        /// TODO: make more efficient by drawing edges directly when construction the new nodes
        for (lemon::ListDigraph::NodeIt r1(p_region_graph); r1 != lemon::INVALID; ++r1)
        {
               for (lemon::ListDigraph::NodeIt r2(p_region_graph); r2 != lemon::INVALID; ++r2)
                {
                        if ((r1 != r2) && (p_region_layer[r1] == p_region_layer[r2]-1))
                        {
                                gbp::VarSet intersection = p_region_qubits[r1] & p_region_qubits[r2];
                                if (!intersection.empty())
                                {
                                        lemon::ListDigraph::Arc new_edge = p_region_graph.addArc(r1, r2);
                                }
                        }
                        
                }
        }



        // fill belief_constituency_list

        // first get all descendants and construct shadow
        lemon::ListDigraph::NodeMap< xt::xarray<int> > shadow(p_region_graph);

        // construct depth-first search
        lemon::Dfs<lemon::ListDigraph> dfs(p_region_graph);

        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                dfs.run(region);
                std::vector<int> desc;
                for (lemon::ListDigraph::NodeIt r2(p_region_graph); r2 != lemon::INVALID; ++r2)
                {
                        if ((dfs.reached(r2) == true) && (r2 != region))
                        {
                                desc.push_back(p_region_graph.id(r2));
                        }
                }
                desc.push_back(p_region_graph.id(region));
                shadow[region] = xt::adapt(desc);
        }

         // construct blanket
        lemon::ListDigraph::NodeMap< xt::xarray<int> > blanket(p_region_graph);
        
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                std::cout << "region " << p_region_graph.id(region) << " ";
                std::vector<int> bcl;
                std::vector<int> blank;
                for (auto s: shadow[region])
                {
                        lemon::ListDigraph::Node s_region = lemon::ListDigraph::nodeFromId(s);
                        for (lemon::ListDigraph::InArcIt edge_from_parent(p_region_graph,s_region); edge_from_parent != lemon::INVALID; ++edge_from_parent)
                        {
                                if (!xt::isin(p_region_graph.id(p_region_graph.source(edge_from_parent)),shadow[region])())
                                {
                                        blank.push_back(p_region_graph.id(p_region_graph.source(edge_from_parent)));
                                        // p_belief_constituency_list[region].push_back(edge_from_parent);
                                        bcl.push_back(p_region_graph.id(edge_from_parent));
                                }
                        }
                }
                blanket[region] = xt::adapt(blank);
                auto last = std::unique(bcl.begin(), bcl.end());
                bcl.erase(last, bcl.end());
                for (auto b : bcl)
                {
                        std::cout << b << " ";
                        p_belief_constituency_list[region].push_back(lemon::ListDigraph::arcFromId(b));
                }
                std::cout << std::endl;
        }


        // for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        // {
        //         // messages from parents
        //         for (lemon::ListDigraph::InArcIt edge_from_parent(p_region_graph,region); edge_from_parent != lemon::INVALID; ++edge_from_parent)
        //         {
        //                 p_belief_constituency_list[region].push_back(edge_from_parent);
        //         }
        //          // messages into children (surface, xzzx: no other descendants) from other parents
        //         for (lemon::ListDigraph::OutArcIt edge_to_child(p_region_graph,region); edge_to_child != lemon::INVALID; ++edge_to_child)
        //         {
        //                 for (lemon::ListDigraph::InArcIt edge_from_parent_of_child(p_region_graph,p_region_graph.target(edge_to_child)); edge_from_parent_of_child != lemon::INVALID; ++edge_from_parent_of_child)
        //                 {
        //                         if (p_region_graph.id(region) != p_region_graph.id(p_region_graph.source(edge_from_parent_of_child)))
        //                         {
        //                                 p_belief_constituency_list[region].push_back(edge_from_parent_of_child);
        //                         }
        //                 }
        //         }

        // }

        // fill message_constituency_list
        for (lemon::ListDigraph::ArcIt edge(p_region_graph); edge != lemon::INVALID; ++edge)
        {
                std::set<int> ids = {};
                ids.insert(p_region_graph.id(edge));
                for (auto e : p_belief_constituency_list[p_region_graph.target(edge)])
                {
                        auto search = ids.find(p_region_graph.id(e));
                        if (search == ids.end())
                        {
                                ids.insert(p_region_graph.id(e));
                        }
                }
                for (auto e : p_belief_constituency_list[p_region_graph.source(edge)])
                {
                        auto search = ids.find(p_region_graph.id(e));
                        if (search == ids.end())
                        {
                                p_message_constituency_list[edge].emplace_back(e);
                                ids.insert(p_region_graph.id(e));
                        }
                }
        }

        // calculate counting numbers
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                if (p_region_layer[region] == 0)
                {
                        p_region_counting_number[region] = 1.0;
                }
                else
                {
                        p_region_counting_number[region] = -1.0;
                }
        }

        // print graph
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                std::cout << "/ region " << p_region_graph.id(region) << std::endl;
                std::cout << "| region_qubits " << p_region_qubits[region] << std::endl;
                std::cout << "| p_region_layer " << p_region_layer[region] << std::endl;
                std::cout << "| out_edges:     ";
                        for (lemon::ListDigraph::OutArcIt edge_to_child(p_region_graph,region); edge_to_child != lemon::INVALID; ++edge_to_child)
                        {
                                std::cout << "--> " <<  p_region_graph.id(p_region_graph.target(edge_to_child)) << " ;";
                        }
                        std::cout << "\n";
                std::cout << "| p_belief_constituency_list[region] ";
                for (auto e : p_belief_constituency_list[region])
                {
                        std::cout << p_region_graph.id(p_region_graph.source(e)) << " -> " <<  p_region_graph.id(p_region_graph.target(e)) << " ";
                } 
                std::cout << std::endl;

                std::cout << "\\" << std::endl;
        }

}

void gbp::Regiongraph::initialize_data(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome)
{
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                if (p_region_layer[region] > 0)
                {
                        p_region_belief[region] = p_region_factor[region];
                        // for (size_t i = 0; i < 2; i++)
                        // {
                        //         p_region_belief[region].set(i, error_probabilities(i));
                        // }
                        // p_region_localfactors_lnsum[region] = p_region_belief[region].log(true);

                        for (auto q : p_region_qubits[region])
                        {
                                gbp::TFactor<long double> qF(q);
                                for (size_t i = 0; i < 2; i++)
                                {
                                        qF.set(i, error_probabilities(i));
                                }
                                p_region_belief[region] *= qF;
                                p_region_localfactors_lnsum[region] += qF.log(true);
                        }
                        // std::cout << "** p_region_graph.id(region) " << p_region_graph.id(region) << "\n" <<  "p_region_belief[region] = " <<  p_region_belief[region] << std::endl;
                }
                else if (p_region_layer[region] == 0)
                {
                        p_region_belief[region] = p_region_factor[region];
                        p_region_belief[region].fill(1.0);
                        p_region_localfactors_lnsum[region] = p_region_factor[region];
                        p_region_localfactors_lnsum[region].fill(0.0);
                        for (auto q : p_region_qubits[region])
                        {
                                gbp::TFactor<long double> qF(q);
                                for (size_t i = 0; i < 2; i++)
                                {
                                        qF.set(i, error_probabilities(i));
                                }
                                p_region_belief[region] *= qF;
                                p_region_localfactors_lnsum[region] += qF.log(true);
                        }
                        // std::cout << "** p_region_graph.id(region) " << p_region_graph.id(region) << "\n" <<  "p_region_belief[region] = " <<  p_region_belief[region] << std::endl;
                        if (p_cvm)
                        {
                                xt::xarray<int> sydromes = xt::index_view(syndrome,p_region_checks[region]);
                                parity_multicheck(p_region_checks[region],sydromes,p_region_belief[region]);
                                // std::cout << "region  " << p_region_graph.id(region) << std::endl;
                                // std::cout << "p_region_checks[region] = " << p_region_checks[region] << std::endl;
                                // std::cout << "sydromes = " << sydromes << std::endl;
                                // std::cout << "p_region_belief[region] = " << p_region_belief[region] << std::endl;
                        }
                        else
                        {
                                filter_parity(syndrome(p_region_checks[region](0)),p_region_belief[region]);
                        }
                        

                }
                p_region_localfactors_product[region] = p_region_belief[region];

        }

        for (lemon::ListDigraph::ArcIt edge(p_region_graph); edge != lemon::INVALID; ++edge)
        {
                p_message[edge] = p_region_factor[p_region_graph.target(edge)];
                p_message[edge].fill(1.0);
        }
        lemon::mapFill(p_region_graph,p_message_converged,false);

}

void gbp::Regiongraph::update_messages(const long double w, bool normalize)
{
        for (lemon::ListDigraph::ArcIt edge(p_region_graph); edge != lemon::INVALID; ++edge)
        {
                if (!p_message_converged[edge])
                {

                        gbp::TFactor<long double> old_message = p_message[edge];
                        gbp::TFactor<long double> belief_marginalization =  p_region_belief[p_region_graph.source(edge)].marginal(p_region_qubits[p_region_graph.target(edge)]) * p_region_belief[p_region_graph.target(edge)].inverse();// - 1.0;
                        
                        gbp::TFactor<long double> new_message =old_message*belief_marginalization;
                        // gbp::TFactor<long double> new_message =old_message*(belief_marginalization*w + 1.0);
                        // new_message.normalize();

                        // gbp::TFactor<long double> old_message = p_message[edge];
                        // gbp::TFactor<long double> belief_marginalization =  p_region_belief[p_region_graph.source(edge)].marginal(p_region_qubits[p_region_graph.target(edge)]) * p_region_belief[p_region_graph.target(edge)].inverse() ;
                        

                        // gbp::TFactor<long double> new_message =old_message.operator^(1.0-w) * belief_marginalization.operator^(w);
                        //
                        if (normalize)
                                new_message.normalize();
                        p_message[edge] = new_message;

                        // long double diff = new_message.get(0) - old_message.get(0);

                        // if (abs(diff) < 1E-8)
                        // {
                        //         p_message_converged[edge] = true;
                        // }
                        // else
                        // {
                        //         p_message[edge] = new_message;
                        // }
                }
        }
}

void gbp::Regiongraph::update_beliefs(bool normalize)
{
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                // local factors
                p_region_belief[region] = p_region_localfactors_product[region];
                // messages
                for (auto edge : p_belief_constituency_list[region])
                {
                        p_region_belief[region] *= p_message[edge];
                }

                if (normalize)
                        p_region_belief[region].normalize();
        }
}

xt::xarray<long double> gbp::Regiongraph::get_tdqs(bool update)
{
        if (update)
        {
                long double overall_avg_energy = 0;
                long double overall_entropy = 0;
                long double overall_free_energy = 0;

                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        gbp::TFactor<long double> belief_normalized = p_region_belief[region].normalized(); /// TODO: Understand why not-normalized is going down well while normalized isn't
                        // region average energy
                        long double region_avg_energy = -(belief_normalized*p_region_localfactors_lnsum[region]).sum();
                        // region entropy
                        long double region_entropy = belief_normalized.entropy();
                        // std::cout << "region " << p_region_graph.id(region) << " belief_normalized = " << belief_normalized << " entropy = " << region_entropy << "\n";
                        // region free energy
                        long double region_free_energy = region_avg_energy-region_entropy;
                        
                        p_region_tdqs[region](0) = region_avg_energy;
                        overall_avg_energy += (p_region_counting_number[region] * region_avg_energy);
                        p_region_tdqs[region](1) = region_entropy;
                        overall_entropy += (p_region_counting_number[region] * region_entropy);
                        p_region_tdqs[region](2) = region_free_energy;
                        overall_free_energy += (p_region_counting_number[region] * region_free_energy);

                }
                p_tdqs(0) = overall_avg_energy;
                p_tdqs(1) = overall_entropy;
                // std::cout << "overall_entropy = " << overall_entropy << "\n";
                p_tdqs(2) = overall_free_energy;
        }
        return p_tdqs;
}

xt::xarray<int> gbp::Regiongraph::get_hard_decision(int method)
{
        p_hard_decision.fill(0);
        if (method == 1)
        {
                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        if (p_region_layer[region] == 0)
                        {
                                // std::cout << "p_region_belief[region].argMax() = " << p_region_belief[region].argMax() << " \n" ;
                                for (const auto& n : p_region_belief[region].argMax())
                                {
                                        if (n.second == 1)
                                        {
                                                // std::cout << n.first << " \n" ;
                                                p_hard_decision(n.first.label()) =  1;
                                        }
                                }
                        }
                }
                p_incompatibility_score = 0;
                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        if (p_region_layer[region] == 0)
                        {
                                int v = p_region_qubits[region].front().label();
                                if (p_hard_decision(v) == 1)
                                {
                                        p_incompatibility_score++;
                                }
                                else if (p_hard_decision(v) >= 2)
                                {
                                        p_hard_decision(v) = 1;
                                }
                        }
                }
        }
        else if (method == 0)
        {
                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        if (p_region_layer[region] == 1)
                        {
                                for (const auto& n : p_region_belief[region].argMax())
                                {
                                        if (n.second == 1)
                                        {
                                                // std::cout << n.first << " \n" ;
                                                p_hard_decision(n.first.label()) =  1;
                                        }
                                }
                        }
                }
        }
        // for (size_t v = 0; v < p_n_qubits; v++)
        // {
        //         if (p_hard_decision(v) ==  1)
        //         {
        //                 p_incompatibility_score++;
        //         }
        //         else if (p_hard_decision(v) ==  2)
        //         {
        //                 p_hard_decision(v) = 1;
        //         }
        // }
        return p_hard_decision;
}

void gbp::Regiongraph::filter_parity(int parity, gbp::TFactor<long double> &belief, bool normalize)
{
        for (size_t linearState = 0; linearState < belief.nrStates(); linearState++)
        { // for all (joint) states
                // calculate states of x0 and x1 corresponding to state linearState of X
                std::map<gbp::Var, size_t> states = gbp::calcState(belief.vars(), linearState);

                const std::size_t state_parity = std::accumulate(std::begin(states), std::end(states), 0, [](const std::size_t previous,  const auto& element) { return previous + element.second; }) % 2;

                if (state_parity != parity)
                {
                        belief.set(linearState,0.0);
                }
        }
        if (normalize)
        {
                belief.normalize();
        }
}

void gbp::Regiongraph::construct_check_factors()
{
        for (size_t c = 0; c < p_n_checks; c++)
        {
                std::vector<gbp::Var> check_vars;
                for (size_t q = 0; q < p_n_qubits; q++)
                {
                        if (p_H(c, q) != 0)
                        {
                                check_vars.push_back(gbp::Var(q, 2));
                        }
                }
                gbp::VarSet check_qubits;
                check_qubits = gbp::VarSet(check_vars.begin(), check_vars.end(), check_vars.size());

                gbp::TFactor<long double> even(check_qubits,1);
                gbp::TFactor<long double> odd(check_qubits,1);

                filter_parity(0,even,false);
                filter_parity(1,odd,false);

                std::vector<gbp::TFactor<long double>> check_parity_vec{even,odd};
                p_check_factors.push_back(check_parity_vec);
        }
}

void gbp::Regiongraph::parity_multicheck(xt::xarray<int> checks, xt::xarray<int> parities, gbp::TFactor<long double> &belief, bool normalize)
{
        for (int c = 0; c < checks.size(); c++)
        {
                // std::cout << "belief = " << belief << std::endl;
                // std::cout << "p_check_factors[checks[c]][parities[c]] = " << p_check_factors[checks[c]][parities[c]] << std::endl;

                belief *= p_check_factors[checks[c]][parities[c]];
        }

        if (normalize)
        {
                belief.normalize();
        }
}

void gbp::Regiongraph::save_regiongraph(std::string suffix)
{
        // copy to undirected Graph for planar embedding
        lemon::ListGraph undirected_graph;
        lemon::Undirector<lemon::ListDigraph> ud(p_region_graph);
        lemon::GraphCopy<lemon::Undirector<lemon::ListDigraph>, lemon::ListGraph> graphCopy(ud, undirected_graph);
        graphCopy.run();

        typedef lemon::dim2::Point<int> Point;
        lemon::ListGraph::NodeMap<Point> coords(undirected_graph);
        lemon::ListGraph::NodeMap<int> shapes(undirected_graph, 1);
        lemon::ListGraph::EdgeMap<double> widths(undirected_graph, 0.5);

        lemon::ListGraph::NodeMap<std::string> node_texts(undirected_graph);

        lemon::PlanarEmbedding<const lemon::ListGraph> pe(undirected_graph);
        pe.run();

        lemon::PlanarDrawing<const lemon::ListGraph> pd(undirected_graph);
        pd.run(pe.embeddingMap());

        std::string name = "region_graph_" + suffix + ".eps";

        lemon::graphToEps(undirected_graph, name)
            .coords(pd.coords())
            .edgeWidths(widths)
            .nodeTexts(node_texts)
            .nodeTextSize(1)
            .nodeShapes(shapes)
            .border(200, 10)
            .drawArrows(true)
            .run();
}

void gbp::Regiongraph::print_regiongraph()
{
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                std::cout << "** region " << p_region_graph.id(region) << ":\n";
                std::cout << "   region_checks:  " << p_region_checks[region] << "\n";
                std::cout << "   region_qubits: " << p_region_qubits[region] << "\n";
                
                if (p_region_layer[region] == 0)
                {
                        std::cout << "   out_edges:     ";
                        for (lemon::ListDigraph::OutArcIt edge_to_child(p_region_graph,region); edge_to_child != lemon::INVALID; ++edge_to_child)
                        {
                                std::cout << "--> " << p_region_qubits[p_region_graph.target(edge_to_child)] << " ;";
                        }
                        std::cout << "\n";
                        
                       
                }
                else
                {
                        std::cout << "   in_edges:      ";
                        for (lemon::ListDigraph::InArcIt edge_from_parent(p_region_graph,region); edge_from_parent != lemon::INVALID; ++edge_from_parent)
                        {
                                std::cout << " " << p_region_checks[p_region_graph.source(edge_from_parent)] << " --> ;";
                        }
                        std::cout << "\n";
                }
        }
}

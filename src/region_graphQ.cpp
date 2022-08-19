#include <gbp/region_graphQ.hpp>

// Quaternary implementation

void gbp::QRegiongraph::construct_regiongraph()
{
        // superregions: for every check one region with a factor dependent on all check vars
        for (size_t c = 0; c < p_n_checks; c++)
        {
                std::vector<gbp::Var> check_vars;
                size_t check_degree = xt::sum(xt::row(p_H, c))();
                check_vars.reserve(check_degree);
                for (size_t q = 0; q < p_n_qubits; q++)
                {
                        if (p_H(c, q) != 0)
                        {
                                check_vars.push_back(gbp::Var(q, 4));
                                p_check_types[c] = p_H(c, q);
                        }
                }
                lemon::ListDigraph::Node new_region = p_region_graph.addNode();
                p_region_check[new_region] = c;
                p_region_qubits[new_region] = gbp::VarSet(check_vars.begin(), check_vars.end(), check_vars.size());
                p_region_factor[new_region] = Factor(p_region_qubits[new_region]);
                p_region_layer[new_region] = 0;
        }

        // small regions: all qubits
        for (size_t q = 0; q < p_n_qubits; q++)
        {
                lemon::ListDigraph::Node new_region = p_region_graph.addNode();
                p_region_check[new_region] = -1;
                gbp::VarSet qubit = gbp::VarSet(gbp::Var(q, 4));
                p_region_qubits[new_region] = qubit;
                p_region_factor[new_region] = Factor(p_region_qubits[new_region]);
                p_region_layer[new_region] = 1;
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

                // calculate counting numbers and build parity filter Factors
                if (p_region_layer[region] == 0)
                {
                        p_region_counting_number[region] = 1.0;
                        
                        // p_region_parity_filters[region].reserve(2);
                        // std::cout << "bf" << p_region_parity_filters[region][0] << std::endl;
                        p_region_parity_filters[region].push_back(p_region_factor[region]);
                        p_region_parity_filters[region].push_back(p_region_factor[region]);
                        p_region_parity_filters[region].shrink_to_fit();
                        // std::cout << "af" << p_region_parity_filters[region][0] << std::endl;
                        for (size_t linearState = 0; linearState < p_region_parity_filters[region][0].nrStates(); linearState++)
                        { // for all (joint) states
                                // calculate states of x0 and x1 corresponding to state linearState of X
                                std::map<gbp::Var, size_t> states = gbp::calcState(p_region_parity_filters[region][0].vars(), linearState);
                                
                                int state_parity = 0;
                                std::size_t check_type = p_check_types(p_region_check[region]);
                                for (auto s = std::begin(states); s != std::end(states);++s)
                                {
                                        state_parity += gf4_trace(gf4_mul(check_type,gf4_conj(s->second)));
                                }
                                // const int state_parity = std::accumulate(std::begin(states), std::end(states), 0, [&check_type](const int previous,  const auto& element) { return gf4_trace(gf4_mul(check_type,gf4_conj(previous))^gf4_mul(check_type,gf4_conj(element.second))); }) ;
                                state_parity %= 2;
                                // std::cout << " p = " << parity << " c = " << check_type  << " states = " << states << " state_parity = " << state_parity << "\n";


                                if (state_parity == 1)
                                {
                                        p_region_parity_filters[region][0].set(linearState,0.0);
                                }
                                else if (state_parity == 0)
                                {
                                        p_region_parity_filters[region][1].set(linearState,0.0);
                                }
                        }
                }
                else
                {
                        p_region_counting_number[region] = 1.0;
                        for (lemon::ListDigraph::InArcIt edge_from_parent(p_region_graph,region); edge_from_parent != lemon::INVALID; ++edge_from_parent)
                        {
                               p_region_counting_number[region] -= 1.0;
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

        // // calculate counting numbers
        // for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        // {
        //         if (p_region_layer[region] == 0)
        //         {
        //                 p_region_counting_number[region] = 1.0;
        //         }
        //         else
        //         {
        //                 p_region_counting_number[region] = 1.0;
        //                 for (lemon::ListDigraph::InArcIt edge_from_parent(p_region_graph,region); edge_from_parent != lemon::INVALID; ++edge_from_parent)
        //                 {
        //                        p_region_counting_number[region] -= 1.0;
        //                 }
        //         }
        // }

        

}

void gbp::QRegiongraph::initialize_data(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome)
{
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                if (p_region_layer[region] == 1)
                {
                        p_region_belief[region] = p_region_factor[region];
                        for (size_t i = 0; i < 4; i++)
                        {
                                p_region_belief[region].set(i, error_probabilities(i));
                        }
                        p_region_localfactors_lnsum[region] = p_region_belief[region].log(true);
                }
                else if (p_region_layer[region] == 0)
                {
                        p_region_belief[region] = p_region_factor[region];
                        p_region_belief[region].fill(1.0);
                        p_region_localfactors_lnsum[region] = p_region_factor[region];
                        p_region_localfactors_lnsum[region].fill(0.0);
                        for (auto q : p_region_qubits[region])
                        {
                                Factor qF(q);
                                for (size_t i = 0; i < 4; i++)
                                {
                                        qF.set(i, error_probabilities(i));
                                }
                                p_region_belief[region] *= qF;
                                p_region_localfactors_lnsum[region] += qF.log(true);
                        }
                        std::size_t check_type = p_check_types(p_region_check[region]);
                        // filter_parity(syndrome(p_region_check[region]), check_type, p_region_belief[region]);
                        filter_parity(syndrome(p_region_check[region]),region);

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

void gbp::QRegiongraph::update_messages(const long double damping, bool normalize)
{
        for (lemon::ListDigraph::ArcIt edge(p_region_graph); edge != lemon::INVALID; ++edge)
        {
                // if (!p_message_converged[edge])
                // {

                        Factor old_message = p_message[edge];
                        Factor belief_marginalization =  p_region_belief[p_region_graph.source(edge)].marginal(p_region_qubits[p_region_graph.target(edge)]) * p_region_belief[p_region_graph.target(edge)].inverse();// - 1.0;

                        // Factor new_message =old_message*(belief_marginalization*damping + 1.0);
                        Factor new_message =old_message*belief_marginalization;
                        
                        if (normalize)
                        {
                                new_message.normalize();
                        }
                        // Factor old_message = p_message[edge];
                        // Factor belief_marginalization =  p_region_belief[p_region_graph.source(edge)].marginal(p_region_qubits[p_region_graph.target(edge)]) * p_region_belief[p_region_graph.target(edge)].inverse() ;
                        

                        // Factor new_message =old_message.operator^(1.0-w) * belief_marginalization.operator^(w);
                        // new_message.normalize();


                        // Factor old_message = p_message[edge];
                        // Factor new_message(old_message);
                        // new_message.fill(1);
                        // for (auto edge : p_message_constituency_list[edge])
                        // {
                        //         new_message *= edge;
                        // }
                        // new_message.marginal(p_region_qubits[p_region_graph.target(edge)])

                        long double diff = dist(new_message,old_message,DISTKL);//new_message.get(0) - old_message.get(0);
                        
                        if (abs(diff) < 1E-8)
                        {
                                p_message_converged[edge] = true;
                        }
                        // std::cout << p_message_converged[edge] << " ";
                        // else
                        // {
                                p_message[edge] = new_message;
                        // }
                // }
        }
        // std::cout << std::endl;
}

void gbp::QRegiongraph::update_beliefs(bool normalize)
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
                {
                        p_region_belief[region].normalize();
                }
        }
}

xt::xarray<long double> gbp::QRegiongraph::get_tdqs(bool update, bool normalize)
{
        if (update)
        {
                long double overall_avg_energy = 0;
                long double overall_entropy = 0;
                long double overall_free_energy = 0;

                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        Factor belief_normalized; 
                        if (normalize)
                                belief_normalized =  p_region_belief[region].normalized(); /// TODO: Understand why not-normalized is going down well while normalized isn't
                        else
                                belief_normalized =  p_region_belief[region];
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

xt::xarray<int> gbp::QRegiongraph::get_hard_decision(int method)
{
        p_hard_decision.fill(0);
        p_incompatibility_score = 0;

        if (method == 2)
        {
                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        if (p_region_layer[region] == 0)
                        {
                                for (const auto& n : p_region_belief[region].argMax())
                                {
                                        if (n.second != 0)
                                        {
                                                // p_hard_decision(n.first.label()) =  n.second;
                                                // max_probs(n.first.label()) = max_prob;
                                                
                                                if (p_hard_decision(n.first.label()) != n.second)
                                                {
                                                        p_incompatibility_score++;
                                                }
                                                p_hard_decision(n.first.label()) ^=  n.second;
                                        }
                                }
                        }
                }
        }
        else if (method == 0)
        {
                 xt::xarray<long double> max_probs = xt::zeros<long double>({p_n_qubits});

                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        if (p_region_layer[region] == 1)
                        {
                                // std::cout << p_region_belief[region].argMax() << " \n" ;
                                // long double max_prob = p_region_belief[region].max();
                                for (const auto& n : p_region_belief[region].argMax())
                                {
                                        if (n.second != 0)
                                        {
                                                // if (max_prob >= max_probs(n.first.label()))
                                                // {
                                                        p_hard_decision(n.first.label()) =  n.second;
                                                        // max_probs(n.first.label()) = max_prob;
                                                // }
                                                // p_hard_decision(n.first.label()) ^=  n.second;
                                        }
                                }
                        }
                }
        }
        else if (method == 3)
        {
                xt::xarray<long double> max_probs = xt::zeros<long double>({p_n_qubits});
                
                xt::xarray<bool> visited_qubit = xt::zeros<bool>({p_n_qubits});

                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        if (p_region_layer[region] == 1)
                        {
                                

                                // std::cout << p_region_belief[region].argMax() << " \n" ;
                                // long double max_prob = p_region_belief[region].max();


                                for (const auto& n : p_region_belief[region].argMax())
                                {
                                        //make sure that when visited, only a compatible config is chosen
                                        // if (visited_qubit[n.first.label()])
                                        // {

                                        // }
                                        if (n.second != 0)
                                        {
                                                // if (max_prob >= max_probs(n.first.label()))
                                                // {
                                                        p_hard_decision(n.first.label()) =  n.second;
                                                        // max_probs(n.first.label()) = max_prob;
                                                // }
                                                // p_hard_decision(n.first.label()) ^=  n.second;
                                        }
                                }

                                // for (const auto& q : p_region_qubits[region])
                                // {
                                //         visited_qubit[q] = true;
                                // }
                        }
                }
        }
        else if (method == 1)
        {
                 xt::xarray<long double> max_probs = xt::zeros<long double>({p_n_qubits});

                for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
                {
                        if (p_region_layer[region] == 0)
                        {
                                // std::cout << p_region_belief[region].argMax() << " \n" ;
                                long double max_prob = p_region_belief[region].max();
                                for (const auto& n : p_region_belief[region].argMax())
                                {
                                        if (n.second != 0)
                                        {
                                                if (max_prob >= max_probs(n.first.label()))
                                                {
                                                        p_hard_decision(n.first.label()) =  n.second;
                                                        max_probs(n.first.label()) = max_prob;
                                                }
                                                // p_hard_decision(n.first.label()) ^=  n.second;
                                        }
                                }
                        }
                }
        }
       

        return p_hard_decision;
}

void gbp::QRegiongraph::filter_parity(int parity, int check_type, Factor &belief)
{
        for (size_t linearState = 0; linearState < belief.nrStates(); linearState++)
        { // for all (joint) states
                // calculate states of x0 and x1 corresponding to state linearState of X
                std::map<gbp::Var, size_t> states = gbp::calcState(belief.vars(), linearState);
                
                int state_parity = 0;
                for (auto s = std::begin(states); s != std::end(states);++s)
                {
                        state_parity += gf4_trace(gf4_mul(check_type,gf4_conj(s->second)));
                }
                // const int state_parity = std::accumulate(std::begin(states), std::end(states), 0, [&check_type](const int previous,  const auto& element) { return gf4_trace(gf4_mul(check_type,gf4_conj(previous))^gf4_mul(check_type,gf4_conj(element.second))); }) ;
                state_parity %= 2;
                // std::cout << " p = " << parity << " c = " << check_type  << " states = " << states << " state_parity = " << state_parity << "\n";


                if (state_parity != parity)
                {
                        belief.set(linearState,0.0);
                }
        }
        belief.normalize();
}


void gbp::QRegiongraph::filter_parity(int parity, lemon::ListDigraph::Node &region)
{
        p_region_belief[region] *= p_region_parity_filters[region][parity];
        p_region_belief[region].normalize();
}

void gbp::QRegiongraph::save_regiongraph(std::string suffix)
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

void gbp::QRegiongraph::print_regiongraph()
{
        for (lemon::ListDigraph::NodeIt region(p_region_graph); region != lemon::INVALID; ++region)
        {
                std::cout << "** region " << p_region_graph.id(region) << ":\n";
                std::cout << "   region_check:  " << p_region_check[region] << "\n";
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
                                std::cout << " " << p_region_check[p_region_graph.source(edge_from_parent)] << " --> ;";
                        }
                        std::cout << "\n";
                }
        }
}


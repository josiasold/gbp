#include <gbp/tanner_graph.hpp>

void gbp::Tannergraph::construct_tannergraph()
{
        // check nodes
        for (size_t c = 0; c < p_n_checks; c++)
        {
                lemon::ListGraph::Node new_checknode = p_tannergraph.addNode();
                p_node_index[new_checknode] = c;

                std::vector<gbp::Var> check_vars;
                size_t check_degree = xt::sum(xt::row(p_H, c))();
                check_vars.reserve(check_degree);
                for (size_t q = 0; q < p_n_qubits; q++)
                {
                        if (p_H(c, q) != 0)
                        {
                                check_vars.push_back(gbp::Var(q, 4));
                                p_node_type[new_checknode] = p_H(c, q); // TODO: Adapt to non-CSS
                        }
                }
                
                p_node_qubits[new_checknode] = gbp::VarSet(check_vars.begin(), check_vars.end(), check_vars.size());
                p_node_belief[new_checknode] = gbp::TFactor<long double>(p_node_qubits[new_checknode]);
                p_counting_number[new_checknode] = 1;
        }

        // qubit nodes
        for (size_t q = 0; q < p_n_qubits; q++)
        {
                lemon::ListGraph::Node new_qubitnode = p_tannergraph.addNode();
                p_node_type[new_qubitnode] = -1;
                p_node_index[new_qubitnode] = q;
                gbp::VarSet qubit = gbp::VarSet(gbp::Var(q, 4));
                p_node_qubits[new_qubitnode] = qubit;
                p_node_belief[new_qubitnode] = gbp::TFactor<long double>(p_node_qubits[new_qubitnode]);
                // p_region_layer[new_qubitnode] = 1;
                p_counting_number[new_qubitnode] = 1;
        }

        // draw edges
        // TODO: make more efficient by drawing edges directly when construction the new nodes
        int id = 0;
        for (lemon::ListGraph::NodeIt n1(p_tannergraph); n1 != lemon::INVALID; ++n1)
        {
                if (p_node_type[n1] > 0) // checks
                {
                        for (lemon::ListGraph::NodeIt n2(p_tannergraph); n2 != lemon::INVALID; ++n2)
                        {
                                if (p_node_type[n2] == -1) // qubits
                                {
                                        if (p_H(p_node_index[n1],p_node_index[n2]) != 0)
                                        {
                                                lemon::ListGraph::Edge new_edge = p_tannergraph.addEdge(n1, n2);
                                                p_edge_type[new_edge] = p_H(p_node_index[n1],p_node_index[n2]);
                                                // std::cout << "p_node_index[n1] = " << p_node_index[n1]  << " p_node_index[n2]  = " << p_node_index[n2] << std::endl;
                                                p_edge_H_indices[new_edge] =  {p_node_index[n1], p_node_index[n2]};
                                                // std::cout << "p_edge_H_indices[new_edge](0,1) = " << p_edge_H_indices[new_edge](0) << "," << p_edge_H_indices[new_edge](1) << std::endl;
                                                p_edge_id[new_edge] = id;
                                                id++;
                                                p_counting_number[n2] -=1; 
                                        }
                                        
                                }
                        }
                }
        }

        // for (lemon::ListGraph::NodeIt n1(p_tannergraph); n1 != lemon::INVALID; ++n1)
        // {
        //         std::cout << "node " << p_tannergraph.id(n1) << "\n\t "
        //         << "p_node_index = " << p_node_index[n1] << "\n\t "
        //         << "p_node_type = " << p_node_type[n1] << "\n\t "
        //         << "p_counting_number = " << p_counting_number[n1] << "\n\t ";
        // }

}

void gbp::Tannergraph::initialize_data(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome)
{
        p_initial_error_probabilities = error_probabilities;
        p_syndrome = syndrome;
        for (lemon::ListGraph::NodeIt qubit(p_tannergraph); qubit != lemon::INVALID; ++qubit)
        {
                if (p_node_type[qubit] == -1) // qubits
                {
                        // std::cout << "qubit " << p_node_index[qubit] << "\n";
                        for (lemon::ListGraph::IncEdgeIt inEdge(p_tannergraph,qubit); inEdge != lemon::INVALID; ++inEdge)
                        {
                                long double q0 = error_probabilities(0) + error_probabilities(p_edge_type[inEdge]);
                                p_m_qc[inEdge] = 2*q0 - 1.0;
                                // std::cout << "inEdge " << p_edge_id[inEdge] << " - q0 = " << q0 << "\n";
                        }
                        for (size_t i = 0; i < 4; i++)
                        {
                                p_node_belief[qubit].set(i, error_probabilities(i));
                        }
                        p_node_localfactors_lnsum[qubit] = p_node_belief[qubit].log(true);

                }
                else // checks
                {
                        // p_node_belief[qubit].fill(1.0);
                        // p_node_localfactors_lnsum[qubit] = p_node_belief[qubit].log(true);
                        // for (auto q : p_node_qubits[qubit])
                        // {
                        //         gbp::TFactor<long double> qF(q);
                        //         for (size_t i = 0; i < 4; i++)
                        //         {
                        //                 qF.set(i, error_probabilities(i));
                        //         }
                        //         p_node_belief[qubit] *= qF;
                        //         p_node_localfactors_lnsum[qubit] += qF.log(true);
                        // }
                        // filter_parity(syndrome(p_node_index[qubit]),p_node_type[qubit],p_node_belief[qubit]);
                }
                p_node_initial_belief[qubit] = p_node_belief[qubit];
        }
}

void gbp::Tannergraph::update_messages(const long double w)
{
       // horizontal step = check to bit
       for (lemon::ListGraph::NodeIt check(p_tannergraph); check != lemon::INVALID; ++check)
        {
                if (p_node_type[check] > 0) // checks
                {
                        for (lemon::ListGraph::IncEdgeIt outEdge(p_tannergraph, check); outEdge != lemon::INVALID; ++outEdge)
                        {
                                long double out_message;
                                if (p_syndrome(p_node_index[check]) == 1)
                                {
                                        out_message = -1.0;
                                }
                                else
                                {
                                         out_message = 1.0;
                                }
                                for (lemon::ListGraph::IncEdgeIt inEdge(p_tannergraph, check); inEdge != lemon::INVALID; ++inEdge)
                                {
                                        if (p_edge_id[outEdge] != p_edge_id[inEdge])
                                        {
                                                out_message *= p_m_qc[inEdge];
                                        }
                                }
                                // if (abs(p_m_cq[outEdge] - out_message) < 1E-5) {p_edge_converged[outEdge] = true;}
                                p_m_cq[outEdge] = out_message;
                                // update r
                                p_r[outEdge](0) = pow(0.5*(1.0 + out_message),w);
                                p_r[outEdge](1) = pow(0.5*(1.0 - out_message),w);
                                // p_r[outEdge](0) = 0.5*(1.0 + out_message);
                                // p_r[outEdge](1) = 0.5*(1.0 - out_message);
                        }
                }
        }

       // vertical step = bit to check
       for (lemon::ListGraph::NodeIt qubit(p_tannergraph); qubit != lemon::INVALID; ++qubit)
        {
                if (p_node_type[qubit] == -1) // qubits
                {
                        // std::cout << "qubit " << p_node_index[qubit] << "\n";
                        for (lemon::ListGraph::IncEdgeIt outEdge(p_tannergraph, qubit); outEdge != lemon::INVALID; ++outEdge)
                        {
                                xt::xarray<long double> out_message_p = p_initial_error_probabilities;
                                // std::cout << "out_message_p 0: " << out_message_p << "\n";
                                for (lemon::ListGraph::IncEdgeIt inEdge(p_tannergraph, qubit); inEdge != lemon::INVALID; ++inEdge)
                                {
                                        if (p_edge_id[outEdge] != p_edge_id[inEdge])
                                        {
                                                for (int p = 0; p < 4; p++)
                                                {
                                                        int parity = p_parity_map(p_edge_type[inEdge],p);
                                                        out_message_p(p) *= p_r[inEdge](parity);
                                                }
                                        }
                                }
                                // std::cout << "out_message_p 1: " << out_message_p << "\n";
                                long double q0 = (out_message_p(0) + out_message_p(p_edge_type[outEdge])) * pow(p_r[outEdge](0),1.0/w - 1.0);
                                long double q1 = 0;
                                for (int p = 1; p < 4; p++)
                                {
                                        if (p != p_edge_type[outEdge])
                                        {
                                                q1 += out_message_p(p);
                                        }
                                }
                                // q1 /= p_r[outEdge](1);
                                q1 *= pow(p_r[outEdge](1),1.0/w - 1.0);

                                long double norm = q0 + q1;
                                long double message = (q0 - q1)/norm;
                                if (abs(p_m_qc[outEdge] - message) < 1E-5) {p_edge_converged[outEdge] = true;}

                                p_m_qc[outEdge] = message;
                                
                                // std::cout << "m_qc = " << p_m_qc[outEdge] << "\n";
                        }
                }
        }
}

void gbp::Tannergraph::update_beliefs(int harddecision_method)
{
        p_hard_decision.fill(0);
        for (lemon::ListGraph::NodeIt qubit(p_tannergraph); qubit != lemon::INVALID; ++qubit)
        {
                p_node_belief[qubit] = p_node_initial_belief[qubit];
                if (p_node_type[qubit] == -1) // qubits
                {
                        xt::xarray<long double> belief = p_initial_error_probabilities;

                        for (lemon::ListGraph::IncEdgeIt inEdge(p_tannergraph, qubit); inEdge != lemon::INVALID; ++inEdge)
                        {
                                for (int p = 0; p < 4; p++)
                                {
                                        int parity = p_parity_map(p_edge_type[inEdge],p);
                                        belief(p) *= p_r[inEdge](parity);
                                }
                        }
                        for (int p = 0; p < 4; p++)
                        {
                                p_node_belief[qubit].set(p,belief(p));
                        }
                        if (harddecision_method == 0)
                        {
                                for (const auto& n : p_node_belief[qubit].argMax())
                                {
                                                p_hard_decision(p_node_index[qubit]) =  n.second;
                                }
                        }
                }
        }
        
        if (harddecision_method == 1)
        {
                xt::xarray<long double> max_probs = xt::zeros<long double>({p_n_qubits});

                for (lemon::ListGraph::NodeIt check(p_tannergraph); check != lemon::INVALID; ++check)
                {
                        // p_node_belief[check] = p_initial_error_probabilities;
                        if (p_node_type[check] > 0) // checks
                        {
                                for (lemon::ListGraph::IncEdgeIt inEdge(p_tannergraph, check); inEdge != lemon::INVALID; ++inEdge)
                                {
                                        p_node_belief[check] *= p_node_belief[p_tannergraph.v(inEdge)];
                                }
                        
                                long double max_prob = p_node_belief[check].max();
                                for (const auto& n : p_node_belief[check].argMax())
                                {
                                        if (n.second != 0)
                                        {
                                                // std::cout << n.first << " \n" ;
                                                if (max_prob >= max_probs(n.first.label()))
                                                {
                                                        p_hard_decision(n.first.label()) =  n.second;
                                                        max_probs(n.first.label()) = max_prob;
                                                }
                                        }
                                }
                        }
                }
        }
}

xt::xarray<int> gbp::Tannergraph::get_hard_decision()
{
        return p_hard_decision;
}

// TODO: adapt to non-CSS
void gbp::Tannergraph::filter_parity(int parity, int check_type, gbp::TFactor<long double> &belief)
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

                state_parity %= 2;


                if (state_parity != parity)
                {
                        belief.set(linearState,0.0);
                }
        }
        belief.normalize();
}


xt::xarray<long double> gbp::Tannergraph::get_tdqs(bool update)
{
        if (update)
        {
                long double overall_avg_energy = 0;
                long double overall_entropy = 0;
                long double overall_free_energy = 0;
                long double overall_kld = 0;

                for (lemon::ListGraph::NodeIt node(p_tannergraph); node != lemon::INVALID; ++node)
                {
                        gbp::TFactor<long double> belief_normalized = p_node_belief[node].normalized(); 
                        // region average energy
                        long double region_avg_energy = -(belief_normalized*p_node_localfactors_lnsum[node]).sum();
                        // region entropy
                        long double region_entropy = belief_normalized.entropy();
                        // std::cout << "region " << p_region_graph.id(region) << " belief_normalized = " << belief_normalized << " entropy = " << region_entropy << "\n";
                        // region free energy
                        long double region_free_energy = region_avg_energy-region_entropy;
                        
                        p_node_tdqs[node](0) = region_avg_energy;
                        overall_avg_energy += (p_counting_number[node] * region_avg_energy);
                        p_node_tdqs[node](1) = region_entropy;
                        overall_entropy += (p_counting_number[node] * region_entropy);
                        p_node_tdqs[node](2) = region_free_energy;
                        overall_free_energy += (p_counting_number[node] * region_free_energy);

                        if (p_node_type[node] == -1) // qubits
                        {
                                gbp::TFactor<long double> kld_bs = belief_normalized/p_node_initial_belief[node];
                                kld_bs = kld_bs.log(true);
                                kld_bs *= belief_normalized;
                                long double single_qubit_kld = kld_bs.sum();
                                overall_kld += single_qubit_kld;

                        }

                        

                }
                // overall_kld = log(overall_kld);
                std::cout << "overall kld = " << overall_kld << "\n";
                // p_tdqs(0) = overall_avg_energy;
                p_tdqs(0) = overall_kld;
                p_tdqs(1) = overall_entropy;
                // std::cout << "overall_entropy = " << overall_entropy << "\n";
                p_tdqs(2) = overall_free_energy;
        }
        return p_tdqs;
}

void gbp::Tannergraph::get_subH(xt::xarray<int>& H_sub, xt::xarray<int> &checks_sub, xt::xarray<int> &qubits_sub, int pauli_type)
{
        std::vector<int> subC;
        std::vector<int> subQ;
        // std::vector<std::vector<int> > subE;

        for (lemon::ListGraph::EdgeIt edge(p_tannergraph); edge != lemon::INVALID; ++edge)
        {
                if (!(p_edge_converged[edge]) && (p_edge_type[edge]==pauli_type))
                {
                        subC.push_back(p_edge_H_indices[edge](0));
                        subQ.push_back(p_edge_H_indices[edge](1));
                        // std::cout << "p_edge_H_indices[edge] = " << p_edge_H_indices[edge] << std::endl;
                        // std::cout << "p_edge_H_indices[edge](0,1) = " << p_edge_H_indices[edge](0) << "," << p_edge_H_indices[edge](1) << std::endl;
                        // std::cout << "p_r[edge] = " << p_r[edge] << std::endl;
                }
        }

        xt::xarray<int> subCx = xt::unique(xt::adapt(subC));
        xt::xarray<int> subQx = xt::unique(xt::adapt(subQ));

        xt::xarray<int> subH = xt::view(p_H,xt::keep(subCx),xt::keep(subQx));


        std::cout << "subC = " << subCx << ", "<< subCx.size() << std::endl;
        std::cout << "subQ = " << subQx << ", "<< subQx.size() << std::endl;
        std::cout << "subH = " << subH << ", "<< subH.size() << std::endl;
        H_sub = subH;
        checks_sub = subCx;
        qubits_sub = subQx;
}
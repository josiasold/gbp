#include <gbp/decoderS.hpp>

xt::xarray<int> gbp::DecoderS::decode(const xt::xarray<long double> &error_probabilities,
                                      const xt::xarray<int> &syndrome_0,
                                      std::string channel)
{
        xt::xarray<int> syndrome = syndrome_0;

        xt::xarray<int> error_guess = xt::zeros<int>({p_n_qubits});
        xt::xarray<int> total_error_guess = xt::zeros<int>({p_n_qubits});
        xt::xarray<long double> adapted_error_probabilities(error_probabilities);
        p_took_iterations = 0;
        p_took_repetitions = 0;
        bool done_all = false;

        for (p_took_repetitions = 0; p_took_repetitions < properties.max_repetitions; p_took_repetitions++)
        {
                if (done_all)
                {
                        break;
                }

                if (xt::sum(syndrome)() != 0)
                {
                        // Timer timer0;
                        if (p_took_repetitions > 0)
                        {
                                if (properties.repeat_p0_strategy == 0)
                                {
                                        p_tannergraph.initialize_data(error_probabilities, syndrome);
                                }
                                else if (properties.repeat_p0_strategy == 1)
                                {
                                        int new_total_error_guess = xt::sum(error_guess)();
                                        if (new_total_error_guess > 0)
                                        {
                                                long double p = 1 - error_probabilities(0);
                                                long double p_adapted = abs((long double)xt::sum(syndrome_0^syndrome)()/(long double)p_n_checks);
                                                p_to_xar(p_adapted,adapted_error_probabilities,channel);
                                                // adapted_error_probabilities(0) = 1 - p_adapted;
                                                // for (size_t i = 1; i < 4; i++)
                                                // {
                                                //         adapted_error_probabilities(i) = p_adapted / 3.0;
                                                // }
                                        }
                                        // Timer t;
                                        p_tannergraph.initialize_data(adapted_error_probabilities, syndrome);
                                        // std::cout << "initialize_data " << t.elapsed() << std::endl;
                                }
                                else if (properties.repeat_p0_strategy == 2)
                                {
                                        int new_total_error_guess = xt::sum(error_guess)();
                                        if (new_total_error_guess > 0)
                                        {
                                                long double p = 1 - error_probabilities(0);
                                                long double p_adapted = abs(p - (long double)new_total_error_guess / (long double)p_n_qubits);
                                                p_to_xar(p_adapted,adapted_error_probabilities,channel);
                                                // adapted_error_probabilities(0) = 1 - p_adapted;
                                                // for (size_t i = 1; i < 4; i++)
                                                // {
                                                //         adapted_error_probabilities(i) = p_adapted / 3.0;
                                                // }
                                        }
                                        // Timer t;
                                        p_tannergraph.initialize_data(adapted_error_probabilities, syndrome);
                                        // std::cout << "initialize_data " << t.elapsed() << std::endl;
                                }
                                else if (properties.repeat_p0_strategy == 3)
                                {
                                        long double p_adapted = 0.09;
                                        p_to_xar(p_adapted,adapted_error_probabilities,channel);
                                        p_tannergraph.initialize_data(adapted_error_probabilities, syndrome);
                                        // std::cout << "initialize_data " << t.elapsed() << std::endl;
                                }

                                
                                
                        }
                        else
                        {
                                p_tannergraph.initialize_data(error_probabilities, syndrome);
                        }
                        // std::cout << "initialize_data: " << timer0.elapsed() << std::endl;
                }
                // p_tannergraph.initialize_data(error_probabilities,
                // syndrome);

                for (size_t iteration = 0; iteration < properties.max_iterations; iteration++)
                {
                        // Timer timer;
                        p_tannergraph.update_messages(properties.damping);
                        // std::cout << "message: " << timer.elapsed() << std::endl;
                        p_tannergraph.update_beliefs(properties.hard_decision_method);
                        // std::cout << "update_beliefs: " << timer.elapsed() << std::endl;
                        error_guess = p_tannergraph.get_hard_decision();
                        xt::xarray<int> s = gf4_syndrome(error_guess, p_H);

                        if (properties.save_history)
                        {
                                size_t hist_index = (p_took_repetitions * properties.max_iterations) + iteration;
                                xt::row(p_error_guess_history, hist_index) = total_error_guess ^ error_guess;
                                xt::row(p_syndrome_history, hist_index) = s ^ syndrome;
                                xt::xarray<long double> tdqs = {0.0, 0.0, 0.0};

                                tdqs += p_tannergraph.get_tdqs(true);
                                xt::row(p_tdqs_history, hist_index) = tdqs;

                                // p_incompatibility_score_history(hist_index)
                                // =
                                // p_tannergraph_X.get_incompatibility_score()+p_tannergraph_Z.get_incompatibility_score();
                                long double diff_kld = abs(xt::row(p_tdqs_history, hist_index - 1)(0) - tdqs(0));
                                std::cout << "diff_kld = " << diff_kld << "\n";
                                if ((iteration > 2) && (diff_kld < 1E-2))
                                {
                                        p_took_iterations += iteration;
                                        break;
                                }
                        }

                        if (s == syndrome)
                        {
                                p_took_iterations += iteration;
                                // total_error_guess ^=
                                // error_guess;
                                done_all = true;
                                break;
                        }

                } // for (size_t iteration = 0; iteration <
                // properties.max_iterations; iteration++)
                if (!done_all)
                {
                        p_took_iterations += properties.max_iterations - 1;
                }
                total_error_guess ^= error_guess;
                xt::xarray<int> s = gf4_syndrome(total_error_guess, p_H);
                if (s == syndrome_0)
                {
                        done_all = true;
                        break;
                }
                else
                {
                        syndrome = syndrome_0 ^ s;
                }
        } // for (p_took_repetitions = 0; p_took_repetitions <
        // properties.max_repetitions; p_took_repetitions++)

        return total_error_guess;
}

void gbp::DecoderS::setProperties(const gbp::PropertySet &opts)
{
        if (opts.hasKey("max_iterations"))
                properties.max_iterations = opts.getAs<int>("max_iterations");
        else
                properties.max_iterations = 100;
        if (opts.hasKey("max_repetitions"))
                properties.max_repetitions = opts.getAs<int>("max_repetitions");
        else
                properties.max_repetitions = 100;
        if (opts.hasKey("repeat_p0_strategy"))
                properties.repeat_p0_strategy = opts.getAs<size_t>("repeat_p0_strategy");
        else
                properties.repeat_p0_strategy = 0;
        if (opts.hasKey("repeat_stopping_criterion"))
                properties.repeat_stopping_criterion = opts.getAs<size_t>("repeat_stopping_criterion");
        else
                properties.repeat_stopping_criterion = 0;
        if (opts.hasKey("verbose"))
                properties.verbose = opts.getAs<int>("verbose");
        else
                properties.verbose = 0;
        if (opts.hasKey("damping"))
                properties.damping = opts.getAs<long double>("damping");
        else
                properties.damping = 1.0;
        if (opts.hasKey("normalize"))
                properties.normalize = opts.getAs<long double>("normalize");
        else
                properties.normalize = true;
        if (opts.hasKey("hard_decision_method"))
                properties.hard_decision_method = opts.getAs<int>("hard_decision_method");
        else
                properties.hard_decision_method = 1;
        if (opts.hasKey("save_history"))
                properties.save_history = opts.getAs<bool>("save_history");
        else
                properties.save_history = false;

        // p_regionGraph.set_print_detail(properties.verbose);

        if (properties.save_history)
        {
                int length = properties.max_repetitions * (properties.max_iterations + 1);
                p_error_guess_history = -1 * xt::ones<int>({length, p_n_qubits});
                p_syndrome_history = -1 * xt::ones<int>({length, p_n_checks});
                p_tdqs_history = -1 * xt::ones<int>({length, 3});
                p_incompatibility_score_history = xt::zeros<int>({length});
        }
}

gbp::PropertySet gbp::DecoderS::getProperties() const
{
        gbp::PropertySet opts;
        opts.set("max_iterations", properties.max_iterations);
        opts.set("max_repetitions", properties.max_repetitions);
        opts.set("repeat_p0_strategy", properties.repeat_p0_strategy);
        opts.set("repeat_stopping_criterion", properties.repeat_stopping_criterion);
        opts.set("verbose", properties.verbose);
        opts.set("hard_decision_method", properties.hard_decision_method);
        opts.set("damping", properties.damping);
        opts.set("save_history", properties.save_history);
        return opts;
}

std::string gbp::DecoderS::printProperties() const
{
        std::stringstream s(std::stringstream::out);
        s << "[";
        s << "max_iterations=" << properties.max_iterations << ",";
        s << "max_repetitions=" << properties.max_repetitions << ",";
        s << "repeat_p0_strategy=" << properties.repeat_p0_strategy << ",";
        s << "repeat_stopping_criterion=" << properties.repeat_stopping_criterion << ",";
        s << "verbose=" << properties.verbose << ",";
        s << "hard_decision_method=" << properties.hard_decision_method << ",";
        s << "damping=" << properties.damping << ",";
        s << "save_history=" << properties.save_history << "]";
        return s.str();
}

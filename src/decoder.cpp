#include <gbp/decoder.hpp>

xt::xarray<int> gbp::Decoder::decode(const xt::xarray<long double> &error_probabilities, const xt::xarray<int> &syndrome_0)
{

        xt::xarray<int> syndrome = syndrome_0;

        xt::xarray<int> syndrome_x = xt::view(syndrome, xt::range(0, p_n_checks_X));          // X-syndrome = syndrome of X-stablizers
        xt::xarray<int> syndrome_z = xt::view(syndrome, xt::range(p_n_checks_X, p_n_checks)); // Z-syndrome = syndrome of Z-stablizers

        xt::xarray<int> error_guess_X = xt::zeros<int>({p_n_qubits});
        xt::xarray<int> error_guess_Z = xt::zeros<int>({p_n_qubits});
        xt::xarray<int> error_guess = xt::zeros<int>({p_n_qubits});

        xt::xarray<int> total_error_guess = xt::zeros<int>({p_n_qubits});

        // First initialization

        /// TODO: check wether this is good
        xt::xarray<long double> error_probabilities_X = {1 - (error_probabilities(1) + error_probabilities(3)), error_probabilities(1) + error_probabilities(3)};
        xt::xarray<long double> error_probabilities_Z = {1 - (error_probabilities(2) + error_probabilities(3)), error_probabilities(2) + error_probabilities(3)};

        p_regionGraph_X.initialize_data(error_probabilities_X, syndrome_z); // Z-syndrome has info about X-errors
        p_regionGraph_Z.initialize_data(error_probabilities_Z, syndrome_x); // X-syndrome has info about Z-errors

        bool done_all = false;
        bool done_x = false;
        bool done_z = false;

        p_took_repetitions = 0;
        p_took_iterations = 0;

        for (p_took_repetitions = 0; p_took_repetitions < p_max_repetitions; p_took_repetitions++)
        {
                syndrome_z = xt::view(syndrome, xt::range(p_n_checks_X, p_n_checks));
                syndrome_x = xt::view(syndrome, xt::range(0, p_n_checks_X));
                
                if (xt::sum(syndrome_z)() != 0)
                {
                        p_regionGraph_X.initialize_data(error_probabilities_X, syndrome_z);
                }
                else
                {
                        error_guess_X = xt::zeros<int>({p_n_qubits});
                        done_x = true;
                }
                if (xt::sum(syndrome_x)() != 0)
                {
                        p_regionGraph_Z.initialize_data(error_probabilities_Z, syndrome_x);
                }
                else
                {
                        error_guess_Z = xt::zeros<int>({p_n_qubits});
                        done_z = true;
                }

                // p_regionGraph_X.get_tdqs();
                // p_regionGraph_Z.get_tdqs();

                int same_count = 0;

                for (size_t iteration = 0; iteration < p_max_iterations; iteration++)
                {
                        if (!done_x)
                        {
                                p_regionGraph_X.update_messages(p_w);
                                p_regionGraph_X.update_beliefs();
                                error_guess_X = p_regionGraph_X.get_hard_decision();
                                xt::xarray<int> s_z = gf2_syndrome(error_guess_X, p_H_Z);
                                if (s_z == syndrome_z)
                                {
                                        done_x = true;
                                }
                        }
                        if (!done_z)
                        {
                                p_regionGraph_Z.update_messages(p_w);
                                p_regionGraph_Z.update_beliefs();
                                error_guess_Z = p_regionGraph_Z.get_hard_decision();
                                xt::xarray<int> s_x = gf2_syndrome(error_guess_Z, p_H_X);

                                if (s_x == syndrome_x)
                                {
                                        done_z = true;
                                }
                        }

                        xt::xarray<int> new_error_guess = error_guess_X + 2 * error_guess_Z;
                        if (new_error_guess == error_guess)
                        {
                                same_count++;
                        }
                        else
                        {
                                same_count = 0;
                                error_guess = new_error_guess;
                        }
                        xt::xarray<int> s = gf4_syndrome(error_guess, p_H);
                        if (p_print_detail == 3)
                        {
                                xt::xarray<int> res_s = s^syndrome;
                                std::cout << "* R:" << p_took_repetitions << " I:" << iteration << "\t** " << container_to_string(error_guess,"e_g",true) << " \n\t\t** " << container_to_string(res_s,"r_s",true) << "\n";
                        }
                        if (p_save_history)
                        {
                                size_t hist_index = (p_took_repetitions*p_max_iterations)+iteration;
                                xt::row(p_error_guess_history, hist_index) = total_error_guess^error_guess;
                                xt::row(p_syndrome_history, hist_index) = s^syndrome;
                                xt::xarray<long double> tdqs = {0.0,0.0,0.0};

                                tdqs += p_regionGraph_X.get_tdqs(!done_x);
                                tdqs += p_regionGraph_Z.get_tdqs(!done_z);
                                // tdqs = tdqs + p_regionGraph_Z.get_tdqs(); /// TODO: Check whether this is good
                                xt::row(p_tdqs_history, hist_index) = tdqs;

                                p_incompatibility_score_history(hist_index) = p_regionGraph_X.get_incompatibility_score()+p_regionGraph_Z.get_incompatibility_score();
                        }
                        if (same_count >= 10)
                        {
                                p_took_iterations+=iteration+1;
                                break;
                        }
                        if (s == syndrome)
                        {
                                p_took_iterations+=iteration+1;
                                break;
                        }
                        if (iteration == p_max_iterations-1)
                        {
                                p_took_iterations+=p_max_iterations;
                        }

                } // for (size_t iteration = 0; iteration < p_max_iterations; iteration++)
                total_error_guess  = total_error_guess ^ error_guess;
                xt::xarray<int> s = gf4_syndrome(total_error_guess, p_H);

                if (s == syndrome_0)
                {
                        done_all = true;
                        return total_error_guess;
                }
                else
                {
                        syndrome = syndrome_0 ^ s;
                }
        } // for (p_took_repetitions = 0; p_took_repetitions < p_max_repetitions; p_took_repetitions++)


        return total_error_guess;
}

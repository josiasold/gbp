#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <experimental/filesystem>


#include <valarray>
#include <vector>
#include <set>

#include <gbp/io_tools.hpp>
#include <gbp/la_tools.hpp>
#include <gbp/error_channels.hpp>
#include <gbp/timing.hpp>
#include <gbp/decoder.hpp>

#include "xtensor/xarray.hpp"
#include "xtensor/xnpy.hpp" // load_npy, dump_npy
#include "xtensor/xio.hpp"  // <<
#include "xtensor/xrandom.hpp" // for seeding random generator
#include <xtensor/xindex_view.hpp>

#include "nlohmann/json.hpp"

int main(int argc, char **argv)
{
    // Construct Timer
    Timer timer;

    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " PATH_TO_INPUT_FILE PATH_TO_OUTPUT_DIR" << std::endl;
        return 1;
    }
    // handle input
    std::string PATH_TO_INPUT_FILE = argv[1];
    std::string PATH_TO_OUTPUT_DIR = argv[2];

    std::ifstream json_input_file(PATH_TO_INPUT_FILE);
    nlohmann::json json_input;
    json_input_file >> json_input;

    std::string pathToCodes = json_input.at("path_to_codes");
    std::string code = json_input.at("code");

    std::string pathToH_X = pathToCodes + code + "_hx.npy";
    std::string pathToH_Z = pathToCodes + code + "_hz.npy";

    long double p = json_input.value("p_error",0.05);
    std::string channel = json_input.value("channel","depolarizing");
    int n_errorsamples = json_input.value("n_errorsamples",1000);

    long double p_initial_strategy = json_input.value("p_initial_strategy",-1.0);
    long double w_gbp = json_input.value("w_gbp",0.9);

    static const bool repetition = json_input.value("repetition",true);
    static const int repetition_strategy = json_input.value("repetition_strategy",0); // 0: converged
    static const int max_repetitions = json_input.value("max_repetitions",10);
    static const int max_iterations = json_input.value("max_iterations",35);

    bool save_raw_data = json_input.value("save_raw_data",false);
    bool save_history = json_input.value("save_history",false);
    int print_detail = json_input.value("print_detail",0); // 0: only results, 1: results and intermediate results, 2: failures, 3: details for each run


    bool return_if_success = json_input.value("return_if_success",true);
    bool only_non_converged = json_input.value("only_non_converged",true);




    // Construct output directory and file

    std::string output_suffix = json_input.value("output_suffix","");
    std::string OUTPUT_DIR;
    if (output_suffix == "")
    {
        OUTPUT_DIR = PATH_TO_OUTPUT_DIR + "/" + timer.contruction_time("dir");
    }
    else if (output_suffix == "tmp")
    {
        OUTPUT_DIR = PATH_TO_OUTPUT_DIR + "/tmp";
        if (std::experimental::filesystem::exists(OUTPUT_DIR))
        {
            std::experimental::filesystem::remove_all(OUTPUT_DIR);
        }
    }
    else
    {
        if (!std::experimental::filesystem::exists(PATH_TO_OUTPUT_DIR))
        {
             std::experimental::filesystem::create_directory(PATH_TO_OUTPUT_DIR);
        }
        OUTPUT_DIR = PATH_TO_OUTPUT_DIR + "/" + timer.contruction_time("dir") + "_" + output_suffix;
    }

    std::experimental::filesystem::create_directory(OUTPUT_DIR);
    std::experimental::filesystem::copy(PATH_TO_INPUT_FILE, OUTPUT_DIR);
    std::ofstream OUTPUT_FILE;
    OUTPUT_FILE.open(OUTPUT_DIR + "/" + code + ".out");

    std::cout << "- start time = " << timer.contruction_time() << "\n";
    std::cout << "Input: " << PATH_TO_INPUT_FILE << std::endl;

    OUTPUT_FILE << "- start time = " << timer.contruction_time() << "\n";
    OUTPUT_FILE << "Input: " << PATH_TO_INPUT_FILE << std::endl;

    // Load X/Z Parity Check Matrices
    std::cout << "pathToH_X: " << pathToH_X << std::endl;
    std::cout << "pathToH_Z: " << pathToH_Z << std::endl;

    xt::xarray<int> H_X = xt::load_npy<int>(pathToH_X);
    xt::xarray<int> H_Z = xt::load_npy<int>(pathToH_Z);


    int n_c_X = H_X.shape()[0];
    int n_c_Z = H_Z.shape()[0];
    int n_q = H_X.shape()[1];
    int n_c = n_c_X + n_c_Z;

    // Construct GF(4) parity check matrix

    xt::xarray<int> H = xt::zeros<int>({n_c, n_q});
    auto H_X_view = xt::view(H, xt::range(0, n_c_X), xt::all());
    auto H_Z_view = xt::view(H, xt::range(n_c_X, n_c), xt::all());

    H_X_view = H_X;
    H_Z_view = 2 * H_Z;

    int rank_H_X = gf2_rank(H_X, n_c_X, n_q);
    int rank_H_Z = gf2_rank(H_Z, n_c_Z, n_q);
    int rank_H = gf4_rank(H, n_c, n_q);

   
    // print details on parity check matrix
    if (print_detail > 0)
    {
        std::cout << "rank(H_X) =\t" << rank_H_X << std::endl;
        std::cout << "rank(H_Z) =\t" << rank_H_Z << std::endl;
        std::cout << "rank(H) =\t" << rank_H << std::endl;
        std::cout << "n_q =\t" << n_q << std::endl;
        std::cout << "n_c_X =\t" << n_c_X << std::endl;
        std::cout << "n_c_Z =\t" << n_c_Z << std::endl;
    }
    

    // error probabilities
    xt::xarray<long double> ps;
    if (p == -1) ps = {0.001,0.0025,0.005,0.0075,0.01,0.02,0.03,0.04};
    else if (p == -2) ps = {0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12};
    else if (p == -3) ps = xt::arange<long double>(0.001,0.5,0.001);
    else if (p == -4) ps = {0.15,0.153,0.156,0.159,0.162,0.165,0.168,0.171};
    else if (p == -5) ps = {7.0/n_q,8.0/n_q};
    else if (p == -6) ps = {0.001,0.003,0.007,0.01,0.02,0.03,0.04,0.05};
    else if (p == -7) ps = {0.001, 0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14};
    else if (p == -8) ps = xt::arange<long double>(0.09,0.101,0.003);
    else if (p == -9) ps = {0.125,0.1275,0.13,0.1325,0.135,0.1375,0.14,0.1425,0.145,0.1475,0.15};
    else if (p == -10) ps = xt::arange<long double>(0.0,0.9,0.05);
    else if (p == -11) ps = xt::arange<long double>(0.0,0.51,0.05);
    else ps = {p};

    std::cout << "p\tch_sc\tsci\tsce\tler\tfail\trepeat\titer\tber" << std::endl;
    OUTPUT_FILE  << "p\tch_sc\tsci\tsce\tler\tfail\trepeat\titer\tber" << std::endl;
    
    // initial word.
    xt::xarray<int> x = xt::zeros<int>({n_q});
    if (print_detail == 3)
    {
        std::cout << container_to_string(x,"x",true) << std::endl;
        OUTPUT_FILE << container_to_string(x,"x",true) << std::endl;
    }

    // construct channel
    NoisyChannel noisyChannel;

    // construct Decoder 
    gbp::Decoder Decoder(H_X, H_Z, max_iterations, max_repetitions, w_gbp, save_history, print_detail);

    
    for (size_t i_p = 0; i_p < ps.size(); i_p++)
    {
        // tracked quantities
        int ch_sc = 0;
        int dec_sci = 0;
        int dec_sce = 0;
        int dec_ler = 0;
        int dec_fail = 0;
        int iterations = 0;
        int repeatsplit  = 0;

        long double p_error = ps(i_p);
        long double p_error_for_decoder;
        if (p_initial_strategy == -1)
        {
            p_error_for_decoder = ps(i_p);
        }
        else
        {
             p_error_for_decoder = p_initial_strategy;
        }
        xt::xarray<long double> p_initial;

        // prepare probability vector for error channels
        if (channel == "x")
        {
            p_initial = {1 - p_error_for_decoder, p_error_for_decoder,0,0};
        }
        else if (channel == "xz")
        {
            long double p_xz = 1-sqrt(1-p_error_for_decoder);
            long double p_x = p_xz*(1-p_xz);
            long double p_y = p_xz*p_xz;
            long double p_z = (1-p_xz)*p_xz;
            p_initial = {1 - p_error_for_decoder, p_x,p_z,p_y};
        }
        else
        {
            p_initial = {1 - p_error_for_decoder, p_error_for_decoder / 3.0, p_error_for_decoder / 3.0, p_error_for_decoder / 3.0};
        }

        // sample error channel
        for (size_t i_e = 0; i_e < n_errorsamples; i_e++)
        {
            xt::xarray<int> y = x;
           if (channel == "depolarizing")
            {
                noisyChannel.send_through_pauli_channel(&y,p_error,0);
            }
            else if (channel == "xz")
            {
                noisyChannel.send_through_pauli_channel(&y,p_error,1);
            }
            else if (channel == "x")
            {
                noisyChannel.send_through_pauli_channel(&y,p_error,2);
            }
            else if (channel == "custom")
            {
                // y(35)=2 ; y(36)=3 ; y(38)=3 ; y(41)=2 ; y(58)=1 ; y(62)=2 ; y(70)=2 ; y(71)=3 ; y(73)=1 ; y(81)=2 ; y(84)=1 ; y(85)=3 ; y(86)=2 ; y(90)=1 ; y(101)=3 ; y(104)=3 ; y(106)=3 ; y(112)=2 ; y(118)=2 ; y(123)=1 ; y(128)=3 ; y(129)=1 ; y(131)=2 ; y(137)=3 ; y(143)=1 ; y(163)=1 ; y(167)=1 ; y(188)=3 ; y(193)=3 ; y(194)=3 ; y(197)=1 ; y(205)=2 ; y(212)=2 ; y(217)=2;
                // y(14) = 3;y(16) = 3;y(18) = 1;y(24) = 3;
                // y(10)=1 ; y(12)=3 ; y(16)=3 ; y(19)=2 ; y(21)=3 ; y(24)=3 ; y(29)=3 ; y(39)=1 ; y(40)=3 ; y(41)=1 ; y(71)=1 ; y(73)=1 ; y(74)=1 ; y(81)=2 ; y(90)=1 ; y(96)=3 ; y(108)=3 ; y(111)=3 ; y(117)=1 ; y(120)=2 ; y(127)=1 ; y(129)=3 ; y(137)=1 ; y(143)=1;
                auto string = xt::view(y,xt::range(81,87));
                string = 1;
            }
            else if (channel == "const_weight")
            {
                int weight = (int)(p_error * n_q);
                noisyChannel.const_weight_error_channel(&y, weight);
            }
            else if (channel == "const_weight_restricted")
            {
                int weight = (int)(p_error * n_q);
                int size_classical = (int) (sqrt(n_q-n_c) + sqrt(n_q+n_c))/2;

                noisyChannel.const_weight_error_channel(&y, weight, size_classical, 2);
            }
            else if (channel == "const_weight_restricted2")
            {
                int weight = (int)(p_error * n_q);
                
                int size_classical = (int) (sqrt(n_q-n_c) + sqrt(n_q+n_c))/2;
                int size_classical_T = (int) (sqrt(n_q+n_c) - sqrt(n_q-n_c))/2;

                int base_qubit =0;//size_classical * size_classical;
                noisyChannel.const_weight_error_channel_T(&y, weight, base_qubit, size_classical, 2);
            }
            else
            {
                std::cerr << "Not a valid channel, choose from {'depolarizing', 'xz', 'x', 'custom', 'const_weight', 'const_weight_restricted'}" << std::endl;
            }
            if (print_detail == 3)
            {
                std::cout << container_to_string(y,"y",true) << std::endl;
                OUTPUT_FILE << container_to_string(y,"y",true) << std::endl;
            }
            if (y == x)
            {
                ch_sc++;
                if (print_detail == 3)
                {
                    std::cout << "++ Channel Success"<< std::endl;
                    OUTPUT_FILE << "++ Channel Success"<< std::endl;
                }
            } // if (y == x)
            else
            {
                xt::xarray<int> s_0 = gf4_syndrome(y, H);
                if (print_detail == 3)
                {
                    std::cout << container_to_string(s_0,"s_0",true) << std::endl;
                    OUTPUT_FILE << container_to_string(s_0,"s_0",true) << std::endl;
                }
                xt::xarray<int> error_guess = Decoder.decode(p_initial,s_0);
                
                xt::xarray<int> s = gf4_syndrome(error_guess, H);
                
                if (save_raw_data)
                {
                    Decoder.dump_history(OUTPUT_DIR);
                }
                
                if (s == s_0)
                {
                    iterations += Decoder.took_iterations();
                    repeatsplit += Decoder.took_repetitions();
                    xt::xarray<int> residual_error = error_guess ^ y;
                    if (residual_error == x)
                    {
                        dec_sci++;
                        if (print_detail == 3)
                        {
                            std::cout << "++ GBP Success (i)"<< std::endl;
                            OUTPUT_FILE << "++ GBP Success (i)"<< std::endl;
                        }
                    }
                    else if (gf4_isEquiv(residual_error, H, n_c, n_q))
                    {
                        dec_sce++;
                        if (print_detail == 3)
                        {
                            std::cout << "++ GBP Success (e)"<< std::endl;
                            OUTPUT_FILE << "++ GBP Success (e)"<< std::endl;
                        }
                    }
                    else
                    {
                        dec_ler++;
                        if (print_detail == 3)
                        {
                            std::cout << "-- GBP Logical Error"<< std::endl;
                            OUTPUT_FILE << "-- GBP Logical Error"<< std::endl;
                        }
                    }

                } // if (s == s_0)
                else
                {
                    dec_fail++;
                    if (print_detail == 3)
                    {
                        std::cout << "-- GBP Failure"<< std::endl;
                        OUTPUT_FILE << "-- GBP Failure"<< std::endl;
                    }
                    if (print_detail == 2)
                    {
                        std::cout << "-- GBP Failure"<< std::endl;
                        OUTPUT_FILE << "-- GBP Failure"<< std::endl;
                        std::cout << container_to_string(y,"| y",true) << std::endl;
                        OUTPUT_FILE << container_to_string(y,"| y",true) << std::endl;
                        std::cout << "__" << std::endl;
                        OUTPUT_FILE << "__" << std::endl;
                    }
                } // else <-- if (s == s_0)
            } // else <-- if (y == x)

            if ((print_detail == 1) && ((n_errorsamples >= 100) && (i_e % (int)(n_errorsamples * 0.1) == 0)))// print status to std::cout every 10% of n_errorsamples
            {
                long double ber = (long double)(dec_ler + dec_fail) / (long double)(i_e+1);
                long double avg_repeatsplit = (long double)(repeatsplit) / (long double)(i_e+1-ch_sc);
                long double avg_iterations = (long double)(iterations) / (long double)(i_e+1-ch_sc);
                
                std::cout << p_error << "\t" << ch_sc << "\t" << dec_sci << "\t" << dec_sce << "\t" << dec_ler << "\t" << dec_fail << "\t" << avg_repeatsplit << "\t" << avg_iterations << "\t" << ber << std::endl;
            }

        } // for (size_t i_e = 0; i_e < n_errorsamples; i_e++)
        long double ber = (long double)(dec_ler + dec_fail) / (long double)n_errorsamples;
        long double avg_repeatsplit = (long double)(repeatsplit) / (long double)(n_errorsamples-ch_sc);
        long double avg_iterations = (long double)(iterations) / (long double)(n_errorsamples-ch_sc);
        
        std::cout << p_error << "\t" << ch_sc << "\t" << dec_sci << "\t" << dec_sce << "\t" << dec_ler << "\t" << dec_fail << "\t" << avg_repeatsplit << "\t" << avg_iterations << "\t" << ber << std::endl;
        OUTPUT_FILE << p_error << "\t" << ch_sc << "\t" << dec_sci << "\t" << dec_sce << "\t" << dec_ler << "\t" << dec_fail << "\t" << avg_repeatsplit << "\t" << avg_iterations << "\t" << ber << std::endl;

    } // for (size_t i_p = 0; i_p < ps.size(); i_p++)

    std::cout << "- stop time = " << timer.current_time() << "\ntotal elapsed time = " << timer.elapsed() << " s" << std::endl;
    OUTPUT_FILE << "- stop time = " << timer.current_time() << "\ntotal elapsed time = " << timer.elapsed() << " s" << std::endl;
    OUTPUT_FILE.close();
}

#include <stdio.h>
#include <iostream>
#include <iomanip>
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
#include <gbp/decoderS.hpp>
#include <gbp/decoderB.hpp>
#include <gbp/properties.hpp>

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
    static const size_t repeat_p0_strategy = json_input.value("repeat_p0_strategy",0);
    static const size_t repeat_stopping_criterion = json_input.value("repeat_stopping_criterion",0);
    
    long double damping = json_input.value("damping",0.9);
    long double normalize = json_input.value("normalize",true);

    static const int repetition_strategy = json_input.value("repetition_strategy",0); // 0: converged
    static const int max_repetitions = json_input.value("max_repetitions",10);
    static const int max_iterations = json_input.value("max_iterations",35);

    bool save_raw_data = json_input.value("save_raw_data",false);
    bool save_history = json_input.value("save_history",false);
    int verbose = json_input.value("verbose",0); // 0: only results, 1: results and intermediate results, 2: failures, 3: details for each run
    int hard_decision_method = json_input.value("hard_decision_method",1);

    bool return_if_success = json_input.value("return_if_success",true);
    bool only_non_converged = json_input.value("only_non_converged",true);

    // construct propertiSet
    gbp::PropertySet properties;
    properties.set("max_iterations",max_iterations);
    properties.set("max_repetitions",max_repetitions);
    properties.set("repeat_p0_strategy",repeat_p0_strategy);
    properties.set("repeat_stopping_criterion",repeat_stopping_criterion);
    properties.set("damping",damping);
    properties.set("normalize",normalize);
    properties.set("verbose",verbose);
    properties.set("hard_decision_method",hard_decision_method);
    properties.set("save_history",save_history);

    gbp::PropertySet propertiesPP;
    propertiesPP.set("max_iterations",20);
    propertiesPP.set("max_repetitions",1);
    propertiesPP.set("repeat_p0_strategy",repeat_p0_strategy);
    propertiesPP.set("repeat_stopping_criterion",repeat_stopping_criterion);
    propertiesPP.set("damping",damping);
    propertiesPP.set("normalize",normalize);
    propertiesPP.set("verbose",verbose);
    propertiesPP.set("hard_decision_method",1);
    propertiesPP.set("save_history",save_history);

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
    if (verbose > 0)
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
    else if (p == -7) ps = {0.001, 0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16};
    else if (p == -8) ps = xt::arange<long double>(0.09,0.101,0.003);
    else if (p == -9) ps = {0.125,0.1275,0.13,0.1325,0.135,0.1375,0.14,0.1425,0.145,0.1475,0.15};
    else if (p == -10) ps = xt::arange<long double>(0.0,0.9,0.05);
    else if (p == -11) ps = xt::arange<long double>(0.0,0.51,0.05);
    else ps = {p};

    std::cout << "p\tch_sc\tsci\tsce\tler\tfail\trep_p\trep_s\titer\tber" << std::endl;
    OUTPUT_FILE  << "p\tch_sc\tsci\tsce\tler\tfail\trep_p\trep_s\titer\tber" << std::endl;
    
    // initial word.
    xt::xarray<int> x = xt::zeros<int>({n_q});
    if (verbose == 3)
    {
        std::cout << container_to_string(x,"x",true) << std::endl;
        OUTPUT_FILE << container_to_string(x,"x",true) << std::endl;
    }

    // construct channel
    NoisyChannel noisyChannel;

    // construct ScalarDecoder 
    gbp::DecoderS ScalarDecoder(H, properties);

    size_t repeat_p0 = 0;
    xt::xarray<long double> repeat_p0s_range = xt::linspace<long double>(0.001, 0.499, 100);
    xt::xarray<long double> repeat_p0s = 0;
    
    for (size_t i_p = 0; i_p < ps.size(); i_p++)
    {
        // tracked quantities
        int ch_sc = 0;
        int dec_sci = 0;
        int dec_sce = 0;
        int dec_ler = 0;
        int dec_fail = 0;
        int repeat_p  = 0;
        int iterations = 0;
        int repeatsplit  = 0;

        long double p_error = ps(i_p);
        long double p_error_for_decoder;

        if (p_initial_strategy == -1)
        {
            p_error_for_decoder = ps(i_p);
        }
        else if (p_initial_strategy < 1)
        {
             p_error_for_decoder = p_initial_strategy;
        }
        else if (p_initial_strategy > 1)
        {
            p_error_for_decoder = ps(i_p);
            repeat_p0 = (size_t)p_initial_strategy;
        }
        xt::xarray<long double> p_initial({4});
        p_to_xar(p_error_for_decoder,p_initial,channel);
        // sample error channel
        for (size_t i_e = 0; i_e < n_errorsamples; i_e++)
        {
            xt::xarray<int> y = x;
            if (channel == "custom")
            {
                // y(4)=1; y(23)=3; y(31)=1; y(36)=2; y(39)=2; y(54)=2; y(63)=1; y(71)=1; y(100)=1; y(110)=1; y(132)=2; y(162)=1; y(190)=1; y(191)=3; y(212)=3; y(226)=1; y(228)=3; y(229)=1; y(247)=2; y(256)=3; y(299)=2; y(307)=1; y(316)=1; y(364)=2; y(395)=1;
                y(0)=2; y(32)=2; y(80)=2;
            }
            else
            {
                noisyChannel.send_through_pauli_channel(&y,p_error,channel);
            }
            if (verbose == 3)
            {
                std::cout << container_to_string(y,"y",true) << std::endl;
                OUTPUT_FILE << container_to_string(y,"y",true) << std::endl;
            }
            if (y == x)
            {
                ch_sc++;
                if (verbose == 3)
                {
                    std::cout << "++ Channel Success"<< std::endl;
                    OUTPUT_FILE << "++ Channel Success"<< std::endl;
                }
            } // if (y == x)
            else
            {
                xt::xarray<int> s_0 = gf4_syndrome(y, H);
                if (verbose == 3)
                {
                    std::cout << container_to_string(s_0,"s_0",true) << std::endl;
                    OUTPUT_FILE << container_to_string(s_0,"s_0",true) << std::endl;
                }
                xt::xarray<int> error_guess;
                xt::xarray<int> s;
                if (repeat_p0 > 0)
                    repeat_p0s = xt::random::choice(repeat_p0s_range,repeat_p0-1);
                for (size_t i_rp0 = 0; i_rp0 <= repeat_p0; i_rp0++)
                {
                    if (i_rp0 > 0)
                    {
                        long double pp = repeat_p0s(i_rp0);
                        p_to_xar(pp,p_initial,channel);
                        // p_initial = {1-pp,pp/3.0,pp/3.0,pp/3.0};
                    }
                    error_guess = ScalarDecoder.decode(p_initial,s_0,channel);
                    
                    s = gf4_syndrome(error_guess, H);
                    
                    if (save_raw_data)
                    {
                        ScalarDecoder.dump_history(OUTPUT_DIR);
                    }
                    if (s == s_0)
                    {
                        repeat_p += i_rp0;
                        break;
                    }
                } // for (size_t i_rp0 = 0; i_rp0 <= repeat_p0; i_rp0++)


                // xt::xarray<int> error_guess = ScalarDecoder.decode(p_initial,s_0);
                
                // xt::xarray<int> s = gf4_syndrome(error_guess, H);
                
                if (save_raw_data)
                {
                    ScalarDecoder.dump_history(OUTPUT_DIR);
                }
                
                if (s == s_0)
                {
                    iterations += ScalarDecoder.took_iterations();
                    repeatsplit += ScalarDecoder.took_repetitions();
                    xt::xarray<int> residual_error = error_guess ^ y;
                    if (residual_error == x)
                    {
                        dec_sci++;
                        if (verbose == 3)
                        {
                            std::cout << "++ GBP Success (i)"<< std::endl;
                            OUTPUT_FILE << "++ GBP Success (i)"<< std::endl;
                        }
                    }
                    else if (gf4_isEquiv(residual_error, H, n_c, n_q))
                    {
                        dec_sce++;
                        if (verbose == 3)
                        {
                            std::cout << "++ GBP Success (e)"<< std::endl;
                            OUTPUT_FILE << "++ GBP Success (e)"<< std::endl;
                        }
                    }
                    else
                    {
                        dec_ler++;
                        if (verbose == 3)
                        {
                            std::cout << "-- GBP Logical Error"<< std::endl;
                            OUTPUT_FILE << "-- GBP Logical Error"<< std::endl;
                        }
                    }

                } // if (s == s_0)
                else
                {
                    // postprocessor
                    
                    // get subMatrices
                    xt::xarray<int> HX_sub;
                    xt::xarray<int> Xchecks_sub;
                    xt::xarray<int> Xqubits_sub;
                    ScalarDecoder.get_subH(HX_sub,Xchecks_sub,Xqubits_sub,1);

                    xt::xarray<int> HZ_sub;
                    xt::xarray<int> Zchecks_sub;
                    xt::xarray<int> Zqubits_sub;
                    ScalarDecoder.get_subH(HZ_sub,Zchecks_sub,Zqubits_sub,2);

                    // construct Decoder 
                    if (HX_sub.size()>0)
                    {
                        gbp::DecoderB XDecoder(HX_sub, propertiesPP, true);
                        xt::xarray<int> s_sub = xt::index_view(s,Xchecks_sub);
                        std::cout << "s_sub = " << s_sub << std::endl;  
                        long double pp = xt::sum(s_sub)()/(long double)s_sub.size();
                        p_to_xar(pp,p_initial,channel);
                        std::cout << "p_initial = " << p_initial << std::endl;  
                        xt::xarray<int> error_guess_sub = XDecoder.decode_separate(p_initial,s_sub);
                        xt::xarray<int> mapped_back = xt::zeros<int>({n_q});
                        
                        for (int i = 0; i < error_guess_sub.size(); i++)
                        {
                            if (error_guess_sub(i) == 1)
                            {
                                mapped_back(Xqubits_sub(i)) = 1;
                            }
                        }
                        std::cout << "error_guess_sub = " << error_guess_sub << std::endl; 
                        std::cout << container_to_string(mapped_back,"| mapped_back",true) << std::endl;
                        error_guess ^= mapped_back;
                    }
                    if (HX_sub.size()>0)
                    {
                        gbp::DecoderB ZDecoder(HZ_sub, properties);
                    }
                    
                    s = gf4_syndrome(error_guess, H);
                    if (s == s_0)
                    {
                        xt::xarray<int> residual_error = error_guess ^ y;
                        if (residual_error == x)
                        {
                            dec_sci++;
                            if (verbose == 3)
                            {
                                std::cout << "++ PP Success (i)"<< std::endl;
                                OUTPUT_FILE << "++ PP Success (i)"<< std::endl;
                            }
                        }
                        else if (gf4_isEquiv(residual_error, H, n_c, n_q))
                        {
                            dec_sce++;
                            if (verbose == 3)
                            {
                                std::cout << "++ PP Success (e)"<< std::endl;
                                OUTPUT_FILE << "++ PP Success (e)"<< std::endl;
                            }
                        }
                        else
                        {
                            dec_ler++;
                            if (verbose == 3)
                            {
                                std::cout << "-- PP Logical Error"<< std::endl;
                                OUTPUT_FILE << "-- PP Logical Error"<< std::endl;
                            }
                        }
                    }
                    else
                    {
                    dec_fail++;
                    if (verbose == 3)
                    {
                        std::cout << "-- GBP Failure"<< std::endl;
                        OUTPUT_FILE << "-- GBP Failure"<< std::endl;
                    }
                    if (verbose == 2)
                    {
                        std::cout << "-- GBP Failure"<< std::endl;
                        OUTPUT_FILE << "-- GBP Failure"<< std::endl;
                        std::cout << container_to_string(y,"| y",true) << std::endl;
                        OUTPUT_FILE << container_to_string(y,"| y",true) << std::endl;
                         std::cout << container_to_string(error_guess,"| error_guess",true) << std::endl;
                        OUTPUT_FILE << container_to_string(error_guess,"| error_guess",true) << std::endl;
                        xt::xarray<int> residual_error = error_guess ^ y;
                        std::cout << container_to_string(residual_error,"| residual_error",true) << std::endl;
                        OUTPUT_FILE << container_to_string(residual_error,"| residual_error",true) << std::endl;
                        std::cout << "--" << std::endl;
                        OUTPUT_FILE << "--" << std::endl;
                    }
                    }
                } // else <-- if (s == s_0)
            } // else <-- if (y == x)

            if ((verbose == 1) && ((n_errorsamples >= 100) && (i_e % (int)(n_errorsamples * 0.1) == 0)))// print status to std::cout every 10% of n_errorsamples
            {
                long double ber = (long double)(dec_ler + dec_fail) / (long double)(i_e+1);
                long double avg_repeat_p = (long double)(repeat_p) / (long double)(i_e+1-ch_sc);
                long double avg_repeatsplit = (long double)(repeatsplit) / (long double)(i_e+1-ch_sc);
                long double avg_iterations = (long double)(iterations) / (long double)(i_e+1-ch_sc);
                
                std::cout << p_error << "\t" << ch_sc << "\t" << dec_sci << "\t" << dec_sce << "\t" << dec_ler << "\t" << dec_fail << "\t" << std::fixed << std::setw(3) << std::setprecision(3) << avg_repeat_p << "\t" << avg_repeatsplit << "\t" << avg_iterations << "\t" << ber << std::endl;
            }

        } // for (size_t i_e = 0; i_e < n_errorsamples; i_e++)
        long double ber = (long double)(dec_ler + dec_fail) / (long double)n_errorsamples;
        long double avg_repeatsplit = (long double)(repeatsplit) / (long double)(n_errorsamples-ch_sc);
        long double avg_iterations = (long double)(iterations) / (long double)(n_errorsamples-ch_sc);
        long double avg_repeat_p = (long double)(repeat_p) / (long double)(n_errorsamples-ch_sc);
        
        std::cout << p_error << "\t" << ch_sc << "\t" << dec_sci << "\t" << dec_sce << "\t" << dec_ler << "\t" << dec_fail << "\t"<< std::fixed << std::setw(3) << std::setprecision(3) << avg_repeat_p << "\t"  << avg_repeatsplit << "\t"<< std::fixed << std::setw(3) << std::setprecision(3) << avg_iterations << "\t"<< std::fixed << std::setw(3) << std::setprecision(3) << ber << std::endl;
        OUTPUT_FILE << p_error << "\t" << ch_sc << "\t" << dec_sci << "\t" << dec_sce << "\t" << dec_ler << "\t" << dec_fail << "\t"<< std::fixed << std::setw(3) << std::setprecision(3) << avg_repeat_p << "\t"  << avg_repeatsplit << "\t"<< std::fixed << std::setw(3) << std::setprecision(3) << avg_iterations << "\t"<< std::fixed << std::setw(3) << std::setprecision(3) << ber << std::endl;

    } // for (size_t i_p = 0; i_p < ps.size(); i_p++)

    std::cout << "- stop time = " << timer.current_time() << "\ntotal elapsed time = " << timer.elapsed() << " s" << std::endl;
    OUTPUT_FILE << "- stop time = " << timer.current_time() << "\ntotal elapsed time = " << timer.elapsed() << " s" << std::endl;
    OUTPUT_FILE.close();
}

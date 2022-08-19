#ifndef IO_TOOLS_HPP_
#define IO_TOOLS_HPP_

#include <valarray>
#include <stdlib.h>
#include <iostream>

#include <string>
#include <sstream>

#include "la_tools.hpp"

template <typename T>
inline void print_container(const T &inp)
{
    for (size_t i = 0; i < inp.size(); i++)
    {
        std::cout << inp[i] << " ";
    }
    std::cout << std::endl;
};

template <typename T>
inline void print_container(const T &inp, const char * description)
{
    std::cout << description << " =\n";
    for (size_t i = 0; i < inp.size(); i++)
    {
        std::cout << inp[i] << " ";
    }
    std::cout << std::endl << std::endl;;
};

template <typename T>
inline void print_container(const T &inp, const char * description, const bool weight)
{
    std::cout << "|" << description << "|" << " = " << hamming_weight(inp) << "\n";
    int c = 0;
    for (size_t i = 0; i < inp.size(); i++)
    {
        if (inp.at(i) != 0)
        {
            std::cout << i << ":" << inp.at(i);
            if (c <  hamming_weight(inp) - 1)
            {
                std::cout << " , ";
            }
            c++;
        }
    }
    std::cout << std::endl;
};


template <typename T>
inline std::string container_to_string(const T &inp, const char * description, const bool weight)
{
    std::stringstream ss;
    int w = hamming_weight(inp);
    ss << description << " |"<<  hamming_weight(inp) <<"|";
    if (w > 0)
    {
        ss << "\t = ";
        int c = 0;
        for (size_t i = 0; i < inp.size(); i++)
        {
            if (inp(i) != 0)
            {
                ss << i << ":" << inp.at(i);
                if (c <  hamming_weight(inp) - 1)
                {
                    ss << " , ";
                }
                c++;
            }
        }
    } 
    
    return ss.str();
};

template <typename T>
inline std::string container_to_string(const T &inp)
{
    std::stringstream ss;
    for (size_t i = 0; i < inp.size(); i++)
    {
        ss << i << ":" << inp.at(i);
        if (i <  inp.size() - 1)
        {
            ss << " , ";
        }
    }

    return ss.str();
};

template <typename T>
inline std::string map_to_string(const std::map<T, size_t> inp)
{
    std::stringstream ss;
    ss << "( ";
    for (const auto& n: inp)
    {
        ss << n.first << ":" << n.second << " ";
    }
    ss << ")";
    return ss.str();
};

// https://stackoverflow.com/questions/15006269/c-get-substring-before-a-certain-char
inline std::string strip_(std::string const& s)
{
    std::string::size_type pos = s.find('_');
    std::string ret;
    if (pos != std::string::npos)
    {
        ret = s.substr(0, pos);
    }
    else
    {
        ret=s;
    }
    return ret;
};

#endif
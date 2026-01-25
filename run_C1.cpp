#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <functional>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <omp.h>
#include "rates.h"



//debugging functions
template<typename T>
void print_vector(const std::vector<T>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i != vec.size() - 1)
            std::cout << ", ";
    }
    std::cout << std::endl;
}

template <typename T>
void print_matrix(const std::vector<std::vector<T>>& matrix) {
    for (const auto& row : matrix) {
        for (const auto& val : row)
            std::cout << val << " ";
        std::cout << '\n';
    }
}


//write to and from csv's
template <typename T>
std::vector<std::vector<T>> write_matrix_from_csv(const std::string &file_name)
{
    std::ifstream file(file_name);
    std::vector<std::vector<T>> result_matrix;
    if (!file.is_open())
    {
        std::cerr << "unable to open file" << std::endl;
        return result_matrix; 
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::vector<T> row;
        std::stringstream lineStream(line);
        std::string cell;

        while (std::getline(lineStream, cell, ','))
        {
            double temp_cell = std::stod(cell);
            row.push_back(static_cast<T>(temp_cell));
        }
        result_matrix.push_back(row);
    }
    file.close();
    return result_matrix;
} 

template <typename T>
void write_csv_from_matrix(const std::string &file_name, const std::vector<std::vector<T>> &data)
{
    std::ofstream outFile(file_name);
     for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            outFile << row[i];
            if (i < row.size() - 1) {
                outFile << ",";
            }
        }
        outFile << "\n";
    }

    outFile.close();
}



//for binomial tau
int find_min_kj(const int N, const std::vector<int> X,
    const std::vector<std::vector<int>> V, int j)
{ 
    std::vector<int> reactions = V[j];
    std::vector<int> possible_vals;
    int min_val;
    for (int i = 0; i < N; ++i)
    {
        if (reactions[i] >= 0)
        {
            continue;
        }
        possible_vals.push_back(std::floor(X[i]/std::abs(reactions[i])));
    }
    if (possible_vals.empty())
    {
        min_val = -1;
    }
    else
    {
        auto min_iterator = std::min_element(possible_vals.begin(), possible_vals.end());
        min_val = *min_iterator;
    }
    
    return min_val;
}

//fixed tau
std::vector<std::vector<int>> binomial_tau(int64_t seed, double tau, double time_step, double end_time, const std::vector<int>& initial_state,
    const std::vector<double>& parameters, 
    std::vector<std::function<double(const std::vector<int>&, const std::vector<double>&)>> rates,
    std::vector<std::vector<int>>& V)
{   
    //initialize seed
    std::mt19937_64 rng(seed);

    //initialize time, mutable state vector
    double time = 0;
    std::vector<int> current_state = initial_state;

    //compute sizes of vectors, initialize current rates vector,
    //binomial vector
    int state_size = initial_state.size();
    int rates_size = rates.size(); 
    std::vector<double> current_rates(rates_size);
    std::vector<int> K(rates_size);

    //results and time record
    std::vector<std::vector<int>> results(static_cast<int>(end_time/time_step) + 1, std::vector<int>(state_size, 0));
    results[0] = initial_state;
    double time_record = time_step;
    int time_index = 1;
    
    while (time < end_time)
    {
        //compute rate functions and set up updated state, step 4
        for (int i = 0; i < rates_size; ++i)
        {
            current_rates[i] = rates[i](current_state, parameters);
        }
        std::vector<int> updated_state = current_state;
        
        //step 5

        time = time + tau;
        if (time > time_record)
        {
            results[time_index] = current_state;
            time_index += 1;
            time_record += time_step;
        }

        
        //steps 6
        for (int j = 0; j < rates_size; ++j)
        {   
            std::vector<int> reaction = V[j];
            int k = find_min_kj(state_size, updated_state, V, j);
            int kj;
            if (k == -1)
            {
                double p = current_rates[j]*tau;
                std::poisson_distribution<int> poiss_dist(p);
                kj = poiss_dist(rng);
            }
            else
            {
                double p = current_rates[j]*tau/k;
                std::binomial_distribution<> binomial_dist(k, p);
                kj = binomial_dist(rng);
            }
           
            K[j] = kj;
            for (int i = 0; i < state_size; ++i)
            {
                if (reaction[i] >= 0)
                {
                    continue;
                }
                updated_state[i] += kj*reaction[i];
            }
        }

        //update state vector, step 7
        for (int j = 0; j < rates_size; ++j)
        {
            std::vector<int> v = V[j];
            for (int i = 0; i < state_size; ++i)
            {
                current_state[i] += v[i]*K[j];
                //std::cout << current_state[i] << std::endl;
            }
        }

    }

    results.back() = current_state;

    return results;

}



int main()
{
    //seed, rng, param set
    int64_t seed = 410;
    std::mt19937_64 rng(seed);
    std::string param_set = "C1.csv";
    std::string param_result = "model_runs_C1.csv";

    //number of runs per parameter
    int runs_per_param = 30;

    //time, tau, time_step
    double end_time = 12;
    double tau = 0.3;
    double time_step = 1;

    //initial state
    int scale = 1;
    std::vector<int> urban_initial_state = {7800*scale, 0, 0, 0, 7800*scale, 2*7800*scale, 0, 0, 2*7800*scale};
    std::vector<int> peri_urban_initial_state = {300*scale, 0, 0, 0, 300*scale, 202*scale, 0, 0, 0, 202*scale, 2*(300 + 202)*scale, 0, 0, 2*(300 + 202)*scale};
    std::vector<int> sylvatic_initial_state = {358*scale, 0, 1*scale, 0, (358 + 1)*scale, 2*(358 + 1)*scale, 0, 0, 2*(358+1)*scale};
    std::vector<int> QOI_initial_state = {0, 0, 0, 0};
    std::vector<int> initial_state;
    initial_state.insert(initial_state.end(), urban_initial_state.begin(), urban_initial_state.end());
    initial_state.insert(initial_state.end(), peri_urban_initial_state.begin(), peri_urban_initial_state.end());
    initial_state.insert(initial_state.end(), sylvatic_initial_state.begin(), sylvatic_initial_state.end());
    initial_state.insert(initial_state.end(), QOI_initial_state.begin(), QOI_initial_state.end());

    //read in state change and parameter matrices
    std::vector<std::vector<int>> V = write_matrix_from_csv<int>("V.csv");
    std::vector<std::vector<double>> params = write_matrix_from_csv<double>(param_set);
    
    //initialize rates
    initialize_rates();

    

    ////initialize results vector
    int num_QoI = 4;
    int num_of_params = params.size();
    std::vector<std::vector<int>> results(num_of_params*runs_per_param, std::vector<int>(num_QoI, 0));


    //loop for model runs

    #pragma omp parallel for
    for (int i = 0; i < num_of_params; ++i)
    {
        std::vector<double> current_param = params[i];
        for (int j = 0; j < runs_per_param; ++j)
        {
            //unique seed for model run
            int64_t current_seed = rng();
            //model run  
            std::vector<std::vector<int>> current_results = binomial_tau(current_seed, tau, time_step, 
                end_time, initial_state, current_param, rates, V);
            //grab end results
            std::vector<int> end_result = current_results.back();
            std::vector<int> QoI(end_result.end() - num_QoI, end_result.end());
            //load to results vec
            results[i*runs_per_param + j] = QoI;
            
        }
    }

    write_csv_from_matrix(param_result, results);



    return 0;
}
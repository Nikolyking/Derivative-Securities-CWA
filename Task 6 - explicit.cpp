#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <iostream>

int main() {
    // -------------------------
    // 1. Parameter Initialization
    // -------------------------
    double S0 = 35.0;     // Initial stock price
    double K = 35.0;      // Strike price
    double T = 1.0;       // Time to maturity (years)
    double r = 0.03;      // Risk-free interest rate
    double q = 0.04;      // Dividend yield
    double sigma = 0.45;  // Volatility

    // -------------------------
    // 2. Grid Setup
    // -------------------------
    double S_max = 3.0 * S0; // Stock price upper bound
    int M = 500;             // Number of stock price steps

    double dS = S_max / M; // Stock price increment
    double dt_stable = std::pow(dS, 2) / (std::pow(sigma, 2) * std::pow(M, 2)); // Stability condition for dt

    // Recompute N based on stability condition
    int N = static_cast<int>(T / dt_stable) + 1;
    double dt = T / N; // Time increment

    std::cout << "Number of time steps (N): " << N << std::endl;

    // -------------------------
    // 3. Stability Check
    // -------------------------
    if (dt > (dS * dS) / (sigma * sigma * M * M)) {
        std::cerr << "Error: Explicit method stability condition violated!" << std::endl;
        return EXIT_FAILURE;
    }

    // -------------------------
    // 4. Grid Initialization
    // -------------------------
    // Initialize stock prices from 0 to S_max
    std::vector<double> stock_prices(M + 1, 0.0);
    for(int i = 0; i <= M; ++i){
        stock_prices[i] = i * dS;
    }

    // Initialize option values grid with zeros
    std::vector<std::vector<double>> option_values(M + 1, std::vector<double>(N + 1, 0.0));

    // Boundary conditions at maturity (t = T)
    for(int i = 0; i <= M; ++i){
        option_values[i][N] = std::max(K - stock_prices[i], 0.0);
    }

    // -------------------------
    // 5. Coefficient Calculation
    // -------------------------
    std::vector<double> alpha(M +1, 0.0);
    std::vector<double> beta(M +1, 0.0);
    std::vector<double> gamma(M +1, 0.0);

    for(int i = 0; i <= M; ++i){
        double j = static_cast<double>(i);
        alpha[i] = 0.5 * dt * (sigma * sigma * j * j - (r - q) * j);
        beta[i] = 1.0 - dt * (sigma * sigma * j * j + (r - q));
        gamma[i] = 0.5 * dt * ((r - q) * j + sigma * sigma * j * j);
    }

    // -------------------------
    // 6. Backward Induction
    // -------------------------
    for(int n = N-1; n >=0; --n){
        // Print progress
        if (n % 100000 == 0){
            std::cout << n << std::endl;
        }
        for(int i = 1; i < M; ++i){
            option_values[i][n] = alpha[i] * option_values[i-1][n+1] + 
                                   beta[i] * option_values[i][n+1] + 
                                   gamma[i] * option_values[i+1][n+1];
            // Early exercise condition for American option
            option_values[i][n] = std::max(option_values[i][n], K - stock_prices[i]);
        }

        // Boundary conditions for stock price extremes
        option_values[0][n] = K;    // When stock price is 0, option value is K
        option_values[M][n] = 0.0; // When stock price is very high, option value is 0
    }

    // -------------------------
    // 7. Extracting the Option Price for S0
    // -------------------------
    // Find the index where stock_prices[i] >= S0
    int initial_price_index = 0;
    for(int i = 0; i <= M; ++i){
        if(stock_prices[i] >= S0){
            initial_price_index = i;
            break;
        }
    }

    double american_put_price = option_values[initial_price_index][0];
    std::cout << "The American put option price at S0 = " << S0 
              << " is approximately: " << std::fixed << std::setprecision(4) 
              << american_put_price << std::endl;

    // -------------------------
    // 8. Data Output for Plotting
    // -------------------------
    std::ofstream outfile("option_values.csv");
    if(!outfile.is_open()){
        std::cerr << "Error: Unable to open file for writing." << std::endl;
        return EXIT_FAILURE;
    }

    // Write CSV header
    outfile << "StockPrice,OptionValue\n";
    for(int i = 0; i <= M; ++i){
        outfile << stock_prices[i] << "," << option_values[i][0] << "\n";
    }
    outfile.close();

    std::cout << "Option values at t=0 have been written to 'option_values.csv' for plotting." << std::endl;

    return EXIT_SUCCESS;
}

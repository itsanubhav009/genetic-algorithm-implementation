#include <iostream>
#include <vector>
#include "GeneticAlgorithm.hpp"

using namespace std;

// Problem definition: weights and values of the items
const vector<int> weights = {2, 3, 4, 5};
const vector<int> values = {3, 4, 5, 6};
const int maxWeight = 5;

// Fitness function for the knapsack problem
double knapsackFitness(const Chromosome& chromosome) {
    int totalWeight = 0;
    int totalValue = 0;

    for (size_t i = 0; i < chromosome.size(); ++i) {
        if (chromosome[i] == 1) {  // If the item is selected
            totalWeight += weights[i];
            totalValue += values[i];
        }
    }

    if (totalWeight > maxWeight) {
        return 0;  // Invalid solution, exceed max weight
    }

    return totalValue;  // Return the total value as the fitness
}

int main() {
    // Genetic Algorithm parameters
    GAParams params = {100, 4, 100, 0.01, 0.7};

    // Create a genetic algorithm instance for the knapsack problem
    GeneticAlgorithm ga(params, knapsackFitness);

    // Run the genetic algorithm
    ga.run();

    // Print the best solution found
    Chromosome bestSolution = ga.getBestSolution();
    double bestFitness = ga.getBestFitness();

    cout << "Best solution found: ";
    for (int gene : bestSolution) {
        cout << gene << " ";
    }
    cout << "\nBest fitness (total value): " << bestFitness << endl;

    return 0;
}

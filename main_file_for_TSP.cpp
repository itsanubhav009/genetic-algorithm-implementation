#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include "GeneticAlgorithm.hpp"

using namespace std;

struct City {
    double x, y;
};

double distance(const City& a, const City& b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

vector<City> cities = {
    {0, 0}, {1, 5}, {5, 2}, {7, 8}, {8, 2},
    {6, 6}, {2, 3}, {4, 7}, {9, 1}, {3, 8}
};

double tspFitness(const Chromosome& chromosome) {
    double totalDistance = 0.0;

    for (size_t i = 0; i < chromosome.size() - 1; ++i) {
        int city1 = chromosome[i];
        int city2 = chromosome[i + 1];
        totalDistance += distance(cities[city1], cities[city2]);
    }

    totalDistance += distance(cities[chromosome.back()], cities[chromosome[0]]);

    return -totalDistance;  // Negative distance because GA maximizes fitness
}

Chromosome orderCrossover(const Chromosome& parent1, const Chromosome& parent2) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, parent1.size() - 1);

    int start = dist(gen);
    int end = dist(gen);

    if (start > end) swap(start, end);

    Chromosome offspring(parent1.size(), -1);

    for (int i = start; i <= end; ++i) {
        offspring[i] = parent1[i];
    }

    int currentIndex = (end + 1) % parent1.size();
    for (int i = 0; i < parent2.size(); ++i) {
        int city = parent2[i];
        if (find(offspring.begin(), offspring.end(), city) == offspring.end()) {
            offspring[currentIndex] = city;
            currentIndex = (currentIndex + 1) % parent1.size();
        }
    }

    return offspring;
}

void swapMutation(Chromosome& chromosome) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, chromosome.size() - 1);

    int i = dist(gen);
    int j = dist(gen);
    swap(chromosome[i], chromosome[j]);
}

class GeneticAlgorithmTSP : public GeneticAlgorithm {
public:
    GeneticAlgorithmTSP(const GAParams& params, std::function<double(const Chromosome&)> fitnessFunction)
        : GeneticAlgorithm(params, fitnessFunction) {}

protected:
    Chromosome crossover(const Chromosome& parent1, const Chromosome& parent2) override {
        return orderCrossover(parent1, parent2);
    }

    void mutate(Chromosome& chromosome) override {
        swapMutation(chromosome);
    }
};

int main() {
    GAParams params = {100, 10, 500, 0.02, 0.8};  // populationSize, chromosomeLength, maxGenerations, mutationRate, crossoverRate

    GeneticAlgorithmTSP ga(params, tspFitness);
    ga.run();

    Chromosome bestSolution = ga.getBestSolution();
    double bestFitness = ga.getBestFitness();

    cout << "Best solution (city order): ";
    for (int city : bestSolution) {
        cout << city << " ";
    }
    cout << "\nBest fitness (negative total distance): " << bestFitness << endl;

    return 0;
}

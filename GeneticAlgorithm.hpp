#ifndef GENETIC_ALGORITHM_HPP
#define GENETIC_ALGORITHM_HPP

#include <vector>
#include <functional>
#include <algorithm>
#include <random>

using Chromosome = std::vector<int>;

struct GAParams {
    int populationSize;
    int chromosomeLength;
    int maxGenerations;
    double mutationRate;
    double crossoverRate;
};

class GeneticAlgorithm {
public:
    GeneticAlgorithm(const GAParams& params, std::function<double(const Chromosome&)> fitnessFunction)
        : params(params), fitnessFunction(fitnessFunction) {
        std::random_device rd;
        gen = std::mt19937(rd());
    }

    void run() {
        initializePopulation();
        for (int generation = 0; generation < params.maxGenerations; ++generation) {
            evaluateFitness();
            select();
            crossover();
            mutate();
        }
        evaluateFitness();  // Final evaluation
    }

    Chromosome getBestSolution() const {
        return bestSolution;
    }

    double getBestFitness() const {
        return bestFitness;
    }

protected:
    virtual Chromosome crossover(const Chromosome& parent1, const Chromosome& parent2) = 0;
    virtual void mutate(Chromosome& chromosome) = 0;

private:
    GAParams params;
    std::function<double(const Chromosome&)> fitnessFunction;
    std::vector<Chromosome> population;
    std::vector<double> fitnessValues;
    std::mt19937 gen;
    Chromosome bestSolution;
    double bestFitness = -std::numeric_limits<double>::infinity();

    void initializePopulation() {
        std::uniform_int_distribution<> dist(0, params.chromosomeLength - 1);
        population.resize(params.populationSize, Chromosome(params.chromosomeLength));
        for (auto& chromo : population) {
            for (int i = 0; i < params.chromosomeLength; ++i) {
                chromo[i] = i;
            }
            std::shuffle(chromo.begin(), chromo.end(), gen);
        }
    }

    void evaluateFitness() {
        bestFitness = -std::numeric_limits<double>::infinity();
        for (const auto& chromo : population) {
            double fitness = fitnessFunction(chromo);
            if (fitness > bestFitness) {
                bestFitness = fitness;
                bestSolution = chromo;
            }
            fitnessValues.push_back(fitness);
        }
    }

    void select() {
        std::vector<Chromosome> newPopulation;
        std::uniform_real_distribution<> dist(0.0, 1.0);

        while (newPopulation.size() < params.populationSize) {
            int idx1 = tournamentSelection();
            int idx2 = tournamentSelection();
            Chromosome child = crossover(population[idx1], population[idx2]);
            newPopulation.push_back(child);
        }
        population = std::move(newPopulation);
    }

    void crossover() {
        // Implement crossover method
    }

    void mutate() {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        for (auto& chromo : population) {
            if (dist(gen) < params.mutationRate) {
                mutate(chromo);
            }
        }
    }

    int tournamentSelection() {
        std::uniform_int_distribution<> dist(0, params.populationSize - 1);
        int bestIdx = dist(gen);
        for (int i = 1; i < 5; ++i) {
            int idx = dist(gen);
            if (fitnessValues[idx] > fitnessValues[bestIdx]) {
                bestIdx = idx;
            }
        }
        return bestIdx;
    }
};

#endif // GENETIC_ALGORITHM_HPP

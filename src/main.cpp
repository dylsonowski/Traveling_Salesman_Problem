#include <array>
#include <utility>
#include <vector>
#include <algorithm>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>

static std::array<std::pair<float, float>, 10> s_cityCords = { {{0, 1}, {2, 3}, {6, 5}, {7, 2.5}, {15, -0.5}, {12, 3.5}, {14, 10}, {9.5, 7.5}, {7.5, 9}, {0.5, 10}} };
//static std::array<std::pair<float, float>, 10> s_cityCords = { {{3, 1}, {2, 4}, {12, 2}, {7, 4.5}, {9, 9}, {3, 1.5}, {16, 11}, {11, 8}, {9, 10}, {2, 7}} };
static unsigned int s_mutationCount = 0;
static std::vector<std::array<float, 3>> experimentConf;
static double s_avgSum = 0;

#define POPULATION_GENE std::array<unsigned int, 10>
int MAX_POPULATION_SIZE; //= { 250, 100, 300, 500 };
float POPULATION_SIZE; //= { 0.8, 0.5, 0.7, 0.9 };
#define NUM_OF_ITERATION 1000
float MUTATION_RATE; //= { 0.2, 0.1, 0.3, 0.5 };

std::vector<std::pair<POPULATION_GENE, double>> CreateKronosAndGaia();
std::pair<POPULATION_GENE, POPULATION_GENE> CreateOffsprings(POPULATION_GENE mother, POPULATION_GENE father);
double CalculateDistance(int startCity, int endCity);
int GenerateRandomNumber(int first, int second);
POPULATION_GENE MutateGene(POPULATION_GENE gene);
void GenerateExperimentsConf();

int main(int argc, char** argv) {	
	GenerateExperimentsConf();
	for (int executeIt = 0; executeIt < 28; executeIt++) {
		if (executeIt == 0) {
			MAX_POPULATION_SIZE = 250;
			POPULATION_SIZE = 0.8 * MAX_POPULATION_SIZE;
			MUTATION_RATE = 0.2;
		}
		else {
			MAX_POPULATION_SIZE = experimentConf.at(executeIt - 1).at(0);
			POPULATION_SIZE = experimentConf.at(executeIt - 1).at(1) * MAX_POPULATION_SIZE;
			MUTATION_RATE = experimentConf.at(executeIt - 1).at(2);
		}

		s_avgSum = 0;
		for (int eventIt = 0; eventIt < 10; eventIt++) {
			std::cout << "Max population size:" << MAX_POPULATION_SIZE << " Population size: " << POPULATION_SIZE << "\n";
			std::vector<std::pair<POPULATION_GENE, double>> population = CreateKronosAndGaia(); //creating parent population
			for (int mainIt = 0; mainIt < NUM_OF_ITERATION; mainIt++) {
				//--- Counting full distance of the gene ---//
				for (auto& distIt : population) {
					double distanceSum = 0;
					for (int geneIt = 0; geneIt < distIt.first.size() - 1; geneIt++) {
						distanceSum += CalculateDistance(distIt.first.at(geneIt) - 1, distIt.first.at(geneIt + 1) - 1);
					}
					distanceSum += CalculateDistance(distIt.first.at(9) - 1, distIt.first.at(0) - 1);
					distIt.second = distanceSum;
				}

				//--- Choosing the best 80% of the parent population ---//
				std::sort(population.begin(), population.end(), [&](std::pair<POPULATION_GENE, double> firstElem, std::pair<POPULATION_GENE, double> secondElem) {
					return firstElem.second < secondElem.second; });
				std::vector<std::pair<POPULATION_GENE, double>> newPopulation;
				int newPopulationSize = 0.8f * population.size();
				s_mutationCount = 0;
				while (newPopulation.size() < newPopulationSize) {
					int motherIt = GenerateRandomNumber(0, newPopulationSize), fatherIt = GenerateRandomNumber(0, newPopulationSize);
					double firstDistance = 0, secondDistance = 0;
					std::pair<POPULATION_GENE, POPULATION_GENE> offsprings = CreateOffsprings(population.at(motherIt).first, population.at(fatherIt).first);

					//--- Applying mutation algorithm ---//
					if (GenerateRandomNumber(0, 100) <= MUTATION_RATE * 100) {
						std::cout << "Mutation!\n";
						s_mutationCount++;
						offsprings.first = MutateGene(offsprings.first);
					}

					if (GenerateRandomNumber(0, 100) <= MUTATION_RATE * 100) {
						std::cout << "Mutation!\n";
						s_mutationCount++;
						offsprings.second = MutateGene(offsprings.second);
					}

					for (int it = 0; it < 9; it++) {
						firstDistance += CalculateDistance(offsprings.first.at(it) - 1, offsprings.first.at(it + 1) - 1);
						secondDistance += CalculateDistance(offsprings.second.at(it) - 1, offsprings.second.at(it + 1) - 1);
					}
					firstDistance += CalculateDistance(offsprings.first.at(9) - 1, offsprings.first.at(0) - 1);
					secondDistance += CalculateDistance(offsprings.second.at(9) - 1, offsprings.second.at(0) - 1);
					newPopulation.emplace_back(std::make_pair(offsprings.first, firstDistance));
					newPopulation.emplace_back(std::make_pair(offsprings.second, secondDistance));
				}

				//--- Filling up population with the best parents ---//
				newPopulation.insert(newPopulation.end(), population.begin(), population.end());
				std::sort(newPopulation.begin(), newPopulation.end(), [&](std::pair<POPULATION_GENE, double> firstElem, std::pair<POPULATION_GENE, double> secondElem) {
					return firstElem.second < secondElem.second; });
				newPopulation.erase(newPopulation.begin() + MAX_POPULATION_SIZE, newPopulation.end());
				population = newPopulation;
				std::cout << "Number of mutation in " << mainIt << " iteration: " << s_mutationCount << "\n";
				system("CLS");
			}

			if (executeIt == 0) {
				std::cout << "Best connection: ";
				for (const auto& conIt : population.at(0).first) {
					std::cout << conIt << " -> ";
				}
				std::cout << population.at(0).first.at(0) << "\nWith cost: " << population.at(0).second << "\nContinue?";
				std::cin.get();
			}
			else {
				s_avgSum += population.at(0).second;
			}

			if (executeIt == 0)
				eventIt = 10;
		}

		if (executeIt != 0) {
			double avgDistance = s_avgSum / 10;
			std::fstream outFile;
			outFile.open("Experiment_Report.csv", std::ios::out | std::ios::app);

			outFile << "Experiment parameters: \nMax Population - " << MAX_POPULATION_SIZE << "\nOffsprings Population size - " << POPULATION_SIZE <<
				"\nMutation rate - " << static_cast<int>(MUTATION_RATE * 100) << "%\nShortest road cost - " << s_avgSum << "\n\n\n";
			outFile.close();
		}
	}
	return 0;
}

std::vector<std::pair<POPULATION_GENE, double>> CreateKronosAndGaia() {
	std::vector<std::pair<POPULATION_GENE, double>> firstPopulation;
	for (int it = 0; it < MAX_POPULATION_SIZE; it++) {
		POPULATION_GENE newGene;
		std::random_device device;
		for (int geneIt = 0; geneIt < newGene.size(); geneIt++) {
			newGene.at(geneIt) = geneIt + 1;
		}
		std::shuffle(newGene.begin(), newGene.end(), std::default_random_engine(device()));
		firstPopulation.emplace_back(std::make_pair(newGene, 0));
	}
	return firstPopulation;
}

std::pair<POPULATION_GENE, POPULATION_GENE> CreateOffsprings(POPULATION_GENE mother, POPULATION_GENE father) {
	try {
		if (mother.size() != father.size())
			throw "Incorrect parents gene size!";

		POPULATION_GENE firstOffspring = father, secondOffspring = mother;
		std::vector<int> usedGenes1Off, usedGenes2Off;
		auto motherIt = std::find(mother.begin(), mother.end(), father.at(0));
		auto fatherIt = std::find(father.begin(), father.end(), mother.at(0));
		for (int it = 0; it < mother.size(); it++) {
			if (std::find(usedGenes1Off.begin(), usedGenes1Off.end(), *motherIt) == usedGenes1Off.end()) {
				firstOffspring.at(std::distance(mother.begin(), motherIt)) = *motherIt;
				usedGenes1Off.push_back(*motherIt);
				motherIt = std::find(mother.begin(), mother.end(), father.at(std::distance(mother.begin(), motherIt)));
			}

			if (std::find(usedGenes2Off.begin(), usedGenes2Off.end(), *fatherIt) == usedGenes2Off.end()) {
				secondOffspring.at(std::distance(father.begin(), fatherIt)) = *fatherIt;
				usedGenes2Off.push_back(*fatherIt);
				fatherIt = std::find(father.begin(), father.end(), mother.at(std::distance(father.begin(), fatherIt)));
			}
		}
		return std::make_pair(firstOffspring, secondOffspring);
	}
	catch (std::string errorMessage) { std::cout << errorMessage << "\n"; }
}

double CalculateDistance(int startCity, int endCity) {
	return sqrt(pow(s_cityCords.at(endCity).first - s_cityCords.at(startCity).first, 2) + pow(s_cityCords.at(endCity).second - s_cityCords.at(startCity).second, 2));
}

int GenerateRandomNumber(int first, int second) {
	std::random_device r_device;
	std::default_random_engine r_engine(r_device());
	std::uniform_int_distribution<int> intDist(first, second);
	return intDist(r_engine);
}

POPULATION_GENE MutateGene(POPULATION_GENE gene) {
	int firstSwap = GenerateRandomNumber(0, gene.size() - 1), secondSwap = GenerateRandomNumber(0, gene.size() - 1);
	unsigned int temp = gene.at(firstSwap);
	gene.at(firstSwap) = gene.at(secondSwap);
	gene.at(secondSwap) = temp;
	return gene;
}

void GenerateExperimentsConf() {
	std::vector<int> maxPopSize = { 100, 300, 500 };
	std::vector<float> popSize = { 0.5f, 0.7f, 0.9f }, mutationRate = { 0.1f, 0.3f, 0.5f };
	for (const auto& maxPopIt : maxPopSize) {
		for (const auto& popSizeIt : popSize) {
			for (const auto& mutationRateIt : mutationRate) {
				experimentConf.push_back({ static_cast<float>(maxPopIt), popSizeIt, mutationRateIt });
			}
		}
	}
}
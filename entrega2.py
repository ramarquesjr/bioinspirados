# coding: utf-8
import random
import math


def func_obj(x):
    n = float(len(x))
    # f_exp = -0.2 * math.sqrt(1/n * sum(np.power(x, 2)))
    t = 0
    for i in range(0, len(x)):
        t += x[i]*x[i]
    f_exp = -0.2 * math.sqrt((1*t)/n)
    t = 0
    for i in range(0, len(x)):
        t += math.cos(2 * math.pi * x[i])
    s_exp = 1/n * t
    f = -20 * math.exp(f_exp) - math.exp(s_exp) + 20 + math.exp(1)
    return f


class GeneticAlgorithm():

    def __init__(self, xmin, xmax, num_gen=2, population=10,
                 cross_rate=1, generations=10, mutate_rate=0.01):
        self.population = []
        self.performance = [0] * population
        self.num_gen = num_gen
        self.xmin = xmin
        self.xmax = xmax
        self.cross_rate = cross_rate
        self.mutate_rate = mutate_rate
        self.population_size = population
        self.generations = generations
        self._create_population(population)

    def _gen_chromo(self):
        chromossome = []
        for x in range(self.num_gen):
            chromossome.append(random.uniform(self.xmin, self.xmax))
        return chromossome

    def _create_population(self, size):
        for i in range(size):
            self.population.append(self._gen_chromo())
        self._evaluate_population()

    def _parse_genes(self, index):
        chromossome = self.population[index]
        answer = []
        for pos in range(self.num_gen):
            answer.append(chromossome[pos])
        return answer

    def _evaluate_population(self):
        for i in range(self.population_size):
            self.performance[i] = self._evaluate_chromossome(i)

    def _evaluate_chromossome(self, index):
        return func_obj(self._parse_genes(index))

    # Métodos de seleção de pais
    def _tournament_sel(self):
        sort = random.sample(range(self.population_size), 2)
        if self.performance[sort[0]] < self.performance[sort[1]]:
            return sort[0]
        else:
            return sort[1]

    def _roulette_sel(self):
        total = 0.0
        for i in self.performance:
            total += i

        roulette = []
        for i in self.performance:
            roulette.append(i/total)

        # Select chromossome index
        parents = []
        for i in range(self.population_size):
            accumulate = 0.0
            x = random.random()
            index = 0
            while accumulate < x:
                accumulate += roulette[index]
                index += 1
            parents.append(index - 1)

        self._blx_alpha_beta(parents)

    def _select(self, method='roulette'):
        if method == 'tournament':
            return self._tournament_sel()
        else:
            return self._roulette_sel()

    # Cruzamento
    def _blx_alpha_beta(self, parents, alpha=0.75, beta=0.25):
        n_pop = []
        index = 0
        while len(n_pop) < self.population_size:
            parent1 = self.population[parents[index]]
            parent2 = self.population[parents[index+1]]
            p1 = [0]*self.num_gen
            p2 = [0]*self.num_gen
            for pos in range(self.num_gen):
                dimension = abs(parent1[pos] - parent2[pos])
                if (dimension < self.xmin and dimension > self.xmax):
                    print(dimension)
                if parent1 <= parent2:
                    print("Entre ", parent1[pos]-alpha*dimension, parent2[pos]+beta*dimension)
                    u = random.uniform(parent1[pos]-alpha*dimension,
                                       parent2[pos]+beta*dimension)
                    p1[pos] = u
                    u = random.uniform(parent1[pos]-alpha*dimension,
                                       parent2[pos]+beta*dimension)
                    p2[pos] = u
                else:
                    print("Entre ", parent2[pos]-beta*dimension, parent1[pos]+alpha*dimension)
                    u = random.uniform(parent2[pos]-beta*dimension,
                                       parent1[pos]+alpha*dimension)
                    p1[pos] = u
                    u = random.uniform(parent2[pos]-beta*dimension,
                                       parent1[pos]+alpha*dimension)
                    p2[pos] = u
            n_pop.append(self._mutate(p1))
            n_pop.append(self._mutate(p2))
            index += 2
        index = self.performance.index(min(self.performance))
        n_pop[0] = self.population[index]
        self.population = n_pop[:]
        self._evaluate_population()

    def _mutate(self, chromo):
        chromossome = chromo
        for pos in range(self.num_gen):
            if random.uniform(0, 1) < self.mutate_rate:
                chromossome[pos] = random.uniform(self.xmin, self.xmax)
        return chromossome

    def view_population(self):
        print('Chromossome numbers: %d' % len(self.population))
        for cont, chromossome in enumerate(self.population):
            print('%2d - %s Fit: %0.2f' % (cont, str(['%.2f' % i for i in chromossome]), self.performance[cont]))

    def best_chromo(self, generation):
        print('\nBest chromossome of %3d generation' % generation)
        index = self.performance.index(min(self.performance))
        chromossome = self.population[index]
        print('%s Fit: %0.2f' % (str(chromossome), self.performance[index]))

    def evolve(self):
        for x in range(self.generations):
            self._select()
            # self.view_population()
            self.best_chromo(x+1)


def main():
    algoritmo = GeneticAlgorithm(-2, 2)
    algoritmo.view_population()
    algoritmo.evolve()
    print("\nPopulação final")
    algoritmo.view_population()


if __name__ == "__main__":
    main()

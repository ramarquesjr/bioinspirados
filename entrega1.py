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


def bin_to_int(binary):
    integer = 0
    for i in range(len(binary)):
        integer += (2**i)*int(binary[i])
    return integer


def mapping(value, xmin, xmax):
    return xmin+((xmax-xmin)/(2**len(value)-1))*bin_to_int(value)


class GeneticAlgorithm():

    def __init__(self, xmin, xmax, num_gen=2, precision=6, population=1000,
                 cross_rate=1, generations=10, mutate_rate=0.01):
        self.population = []
        self.performance = [0] * population
        self.num_gen = num_gen
        self.xmin = xmin
        self.xmax = xmax
        self.cross_rate = cross_rate
        self.mutate_rate = mutate_rate
        self.precision = precision
        self.population_size = population
        self.generations = generations
        self._create_population(population)

    def _gen_chromo(self):
        chromossome = []
        for x in range(self.num_gen*self.precision):
            chromossome.append(random.randint(0, 1))
        return chromossome

    def _create_population(self, size):
        for i in range(size):
            self.population.append(self._gen_chromo())
        self._evaluate_population()

    def _parse_genes(self, index):
        chromossome = self.population[index]
        answer = []
        x = ''
        for pos in range(self.num_gen*self.precision):
            x += str(chromossome[pos])
            if (pos+1) % self.precision == 0:
                answer.append(mapping(x, self.xmin, self.xmax))
                x = ''
        return answer

    def _evaluate_population(self):
        for i in range(self.population_size):
            self.performance[i] = self._evaluate_chromossome(i)

    def _evaluate_chromossome(self, index):
        return func_obj(self._parse_genes(index))

    def _tournament_sel(self):
        sort = random.sample(range(self.population_size), 2)
        if self.performance[sort[0]] < self.performance[sort[1]]:
            return sort[0]
        else:
            return sort[1]

    def _select(self, method='tournament'):
        if method == 'tournament':
            return self._tournament_sel()

    def _cross(self):
        n_pop = []
        while len(n_pop) < self.population_size-1:
            c1 = self._select()
            c2 = self._select()
            while c1 == c2:
                c2 = self._select()
            if random.uniform(0, 1) <= self.cross_rate:
                child1 = self.population[c1][:self.precision] +\
                         self.population[c2][self.precision:self.num_gen*self.precision]
                child2 = self.population[c2][:self.precision] +\
                         self.population[c1][self.precision:self.num_gen*self.precision]
                n_pop.append(self._mutate(child1))
                n_pop.append(self._mutate(child2))
            else:
                n_pop.append(self._mutate(self.population[c1]))
                n_pop.append(self._mutate(self.population[c2]))
        # get best chromossome from population (elitism)
        index = self.performance.index(min(self.performance))
        n_pop[0] = self.population[index]
        self.population = n_pop[:]
        self._evaluate_population()

    def _mutate(self, chromo):
        chromossome = chromo
        for pos in range(self.num_gen*self.precision):
            if random.uniform(0, 1) < self.mutate_rate:
                chromossome[pos] = 1 if chromossome[pos] == 0 else 0
        # print(''.join([str(i) for i in chromo]),
        #       ''.join([str(i) for i in chromossome]))
        return chromossome

    def view_population(self):
        print('Number of chromossome: %d' % len(self.population))
        x = ''
        for cont, chromossome in enumerate(self.population):
            for pos in range(2*self.precision):
                x += str(chromossome[pos])
                if (pos+1) % self.precision == 0:
                    print('%s %4f' % (x, mapping(x, self.xmin, self.xmax)),
                          end=' ')
                    x = ''
            print('   Fitness: %.2f' % self.performance[cont])

    def best_chromo(self, generation):
        print('\nBest chromossome of %3d generation' % generation)
        index = self.performance.index(min(self.performance))
        chromossome = self.population[index]
        x = ''
        for pos in range(2*self.precision):
            x += str(chromossome[pos])
            if (pos+1) % self.precision == 0:
                print('%s %4f' % (x, mapping(x, self.xmin, self.xmax)),
                      end=' ')
                x = ''
        print('   Fitness: %.2f\n' % (self.performance[index]))

    def evolve(self):
        for x in range(self.generations):
            self._cross()
            self.best_chromo(x+1)


def main():
    algoritmo = GeneticAlgorithm(-2, 2)
    algoritmo.evolve()


if __name__ == "__main__":
    main()

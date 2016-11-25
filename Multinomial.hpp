/*

Multinomial distribution

usage

Multinomial mult;
mult.setP(probs,len);           // set vector with cluster probabilities (double
                                   array of length len)
mult.drawSamples(numSamples);   // sample numSamples samples from the
                                   multinomial determined by p (numSamples is an
                                   unsigned long)

Copyright (C) 2016 Maxim Dolgov m<dot>dolgov<at>web<dot>de

Without any warranty.
Free as in speech and free as in beer.

*/

#ifndef MULTINOMIAL_H_
#define MULTINOMIAL_H_

#include <random>
#include "time.h"
#include <stdexcept>

class Multinomial {
 private:
  std::vector<double> probs;
  std::vector<double> cumsum;
  // std::random_device rd;  // produces random integers
  std::mt19937 gen;  // mersenne twister random number generator
  std::uniform_real_distribution<double> dis;  // distribution

 public:
  Multinomial() {
    // gen.seed(rd());  // initialize random number generator unsing a random
    // number from random_device
    gen.seed(time(0));  // initialize random number generator using time
    dis.param(std::uniform_real_distribution<double>::param_type(
        0.0, 1.0));  // set distribution parameters
  };
  ~Multinomial(){};

  void setP(double* p, long len) {
    // copy probabilities
    probs.assign(p, p + len);

    // compute cumulative sum
    cumsum.resize(probs.size());
    // cumulative sum without a positiveness check
    // std::partial_sum(probs.begin(), probs.end(), cumsum.begin());

    // cumulative sum with a positiveness check
    cumsum[0] = probs[0];
    for (int i = 1; i < probs.size(); i++) {
      // check if probabilities are positive and < 1
      if ((probs[i] < 0.0) || (probs[i] > 1.0))
        throw std::out_of_range(
            "Provided vector of probabilities contains negative values or has "
            "entries larger than one.");
      cumsum[i] = cumsum[i - 1] + probs[i];
    }

    // check if the provided vector of probabilities sums to one
    if (abs(cumsum[cumsum.size() - 1] - 1.0) > 1e-8) {
      throw std::out_of_range(
          "Provided vector of probabilities does not sum to one.");
    }
  }

  std::vector<unsigned int> drawSamples(unsigned long numSamples) {
    std::vector<unsigned int> samples;
    samples.reserve(numSamples);

    //#pragma omp parallel for
    for (unsigned long i = 0; i < (unsigned long)numSamples; i++) {
      samples.push_back(drawSample());  // draw a single sample
    }

    return samples;
  }

 private:
  unsigned int drawSample() {
    double randNumber = dis(gen);  // draw sample from the uniform distribution
    unsigned int counter = 0;

    while (counter < (unsigned int)cumsum.size()) {
      if (randNumber < cumsum[counter]) return ++counter;
      ++counter;
    }
    return ++counter;
  }
};

#endif
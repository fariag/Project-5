#include <armadillo>
#include <iostream>
#include "lib.h"
using namespace arma;
using namespace std;

#define   ZERO       1.0E-10

double ran2(long *);

class VMCSolver
	{
	public:
		VMCSolver();
		void runMonteCarloIntegration();
		void setParameters(double alp, double be, double ste);

	private:
		double waveFunction(const mat &r);
		double localEnergy(const mat &r);
		int dimension;
		int charge;
		double stepLength;
		int particles;
		double h;
		double h2;
		long idum;
		double w;
		int mcs;
		mat r;
		double alpha;
		double beta;
		double r12_total;
		mat rOld;
		mat rNew;
	};
	VMCSolver::VMCSolver() :
		dimension(3),
		charge(2),
		particles(2),
		h(0.001),
		h2(1000000),
		idum(-1),
		w(0),
		r12_total(0),
		mcs(1000000)
	{
	}

	void VMCSolver::setParameters(double alp, double be, double om) {
		// Sets parameters for alpha, beta and omega/w
		alpha = alp;
		beta = be;
		w = om;
		if (0.75 <= w && 1 >= w) {
			stepLength = ((9+4.5*alp)/(1+5.8*alp));
		}

		if (0.1 <= w && 0.75 > w) {
			stepLength = ((14+7.5*alp)/(1+6.7*alp));
		}

		if (w < 0.1) {
			stepLength = ((48+7.3*alp)/(1 + 1.7*alp));
		}
		

	}


	void VMCSolver::runMonteCarloIntegration() {
		rOld = zeros<mat>(particles, dimension);
		rNew = zeros<mat>(particles, dimension);
		double wave_old = 0;
		double wave_new = 0;
		double energyTotal = 0;
		double energySquareTotal = 0;
		double delta;
		int accept_count = 0;

		// initial trial positions
		for(int i = 0; i < particles; i++) {
			for(int j = 0; j < dimension; j++) {
				rOld(i,j) = stepLength * (ran2(&idum)-0.5);
			}
		}
		rNew = rOld;
		// loop over Monte Carlo Cycles
		for(int i = 0; i < mcs; i++) {
			
			// Calculate the wave function for the current configuration
			wave_old = waveFunction(rOld);
			
			// Suggest new configuration
			for(int i = 0; i < particles; i++) {
				for(int j = 0; j < dimension; j++) {
					rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum)-0.5);
				}

				// Recalculate the value of the wave function for the new configuration
				wave_new = waveFunction(rNew);
				// Use a random number generator to check if it the new configuration is accepted
				if (ran2(&idum) <= (wave_new*wave_new) / (wave_old*wave_old)) {
					accept_count++;
					for (int j = 0; j < dimension; j++) {
						rOld(i,j) = rNew(i,j);
					}
					wave_old = wave_new;
				} else {
					for(int j = 0; j < dimension; j++) {
						rNew(i,j) = rOld(i,j);
					}
				}
				// update energies
				delta = localEnergy(rNew);
				energyTotal += delta;
				energySquareTotal += delta*delta;
			}
		}
		// Records and prints out results
		double accepted_moves = (double) accept_count/(mcs*particles);
		cout << " " << accepted_moves;
		r12_total = (double) r12_total/(mcs*particles);
		double energy = energyTotal/(mcs * particles);
		double energySquared = energySquareTotal/(mcs * particles);
		double energyVariance = energySquared - energy*energy;
		cout << " " << energy << " " << energyVariance << " " << r12_total << endl;
	}

	double VMCSolver::localEnergy(const mat &r) {
		// Calculates the local energy in the configuration using the analytical expression 

		double localEnergy = 0;
		double rSingleParticle = 0;
		double rr = 0;
		for (int i = 0; i < particles; i++) {
			rSingleParticle = 0;
			for (int j = 0; j < dimension; j++) {
				rSingleParticle += r(i,j)*r(i,j);
			}
			rr += rSingleParticle;
		}

		localEnergy = 0.5*w*w*rr*(1-(alpha*alpha)) + 3*alpha*w;

		double r12 = 0;
		// Finds the electron-electron potential
		for (int i = 0; i < particles; i++) {
			for (int j = i + 1; j < particles; j++) {
				r12 = 0;
				for (int k = 0; k < dimension; k++) {
					r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
				}
				r12 = sqrt(r12);
				r12_total += r12;
				
				localEnergy += (1 / r12);
			}
		}

		return localEnergy;
	}

	double VMCSolver::waveFunction(const mat &r) {
		// Calculates the first wave function for the given configuration
		double argument = 0;
		for (int i = 0; i < particles; i++) {
			double rSingleParticle = 0;
			for (int j = 0; j < dimension; j++) {
				rSingleParticle += r(i,j) * r(i,j);
			}
			argument += rSingleParticle;
		}
		return exp(-argument * alpha * w / 2);
	}

	int main(int argc, char *argv[]) {
		double a_min = stod(argv[1]);
		double a_max = stod(argv[2]);
		double set_w = stod(argv[3]);
		for (double t=a_min; t < a_max+0.01; t += 0.01) {
			printf("%.2f ", t);
			VMCSolver *solver = new VMCSolver();	
			solver->setParameters(t, 0, set_w);
			solver->runMonteCarloIntegration(); 
		}
		return 0;
}

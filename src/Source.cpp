#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#define INPUT_F_NAME "\\input.txt"
#define OUTPUT_F_NAME "\\output.txt"

#define E 2.718281828459

using namespace std;

void findBestLambda(double**, size_t**, double&, double, double, bool(double));
void estimatePhies(double&, const double&, double*, size_t, double);
double estimateResult(size_t*);
size_t* findXs(double**);
bool validate(size_t*);
void input(ifstream&);
void printPhies(double**, ofstream&);
void showInfo();

size_t N, Z, W, C;
double* restrictionsW;
double* restrictionsC;
double* probabilities;

int main() {
	ifstream ifs;
	stringstream ss;
	string wDirectory;
	double** phies = new double* [N];
	size_t* xs = nullptr;
	double lambda;

	showInfo();

	cout << "Please enter the working directory:: ";
	getline(cin, wDirectory);

	ifs.open(wDirectory + string(INPUT_F_NAME));

	if (!ifs.is_open()) {
		return errno;
	}

	input(ifs);
	ifs.close();

	Z = C;
	ss << setprecision(0);

	for (size_t i = 0; i < N; i++) {
		phies[i] = new double[Z];
	}

	// negatives
	findBestLambda(phies, &xs, lambda, -0.002, 0.00001, [](double i) {
		return i < 0.0;
		});

	// positives
	findBestLambda(phies, &xs, lambda, 0.002, -0.00001, [](double i) {
		return i > 0.0;
		});

	ofstream ofs;

	ofs.open(wDirectory + string(OUTPUT_F_NAME));

	if (!ofs.is_open()) {
		return errno;
	}

	ss << "\nLambda = " << lambda << "\nZ | ";
	for (size_t i = 0; i < N; i++) {
		ss << "phi" << i << ", ";
	}
	ss << "\n";

	ofs << ss.str();

	printPhies(phies, ofs);
	ss.str(string());
	ofs.flush();

	ss << setprecision(0);

	ss << "\nXs:: ";
	for (size_t i = 0; i < N; i++) {
		ss << xs[i] << ", ";
	}

	double result = estimateResult(xs);
	double wSum = 0, cSum = 0;

	ss << "\n\nResult:: \nW | ";

	for (size_t i = 0; i < N; i++) {
		wSum += xs[i] * restrictionsW[i];
		ss << xs[i] << "x" << restrictionsW[i] << (i == (N - 1) ? " = " : " + ");
	}
	ss << wSum << "\nC | ";

	for (size_t i = 0; i < N; i++) {
		cSum += xs[i] * restrictionsC[i];
		ss << xs[i] << "x" << restrictionsC[i] << (i == (N - 1) ? " = " : " + ");
	}
	ss << cSum << "\nOverall | " << result;

	ofs << ss.str();
	ofs.close();
	ss.str(string());

	for (size_t i = 0; i < N; i++) {
		delete phies[i];
	}

	delete xs;
	delete[] restrictionsW;
	delete[] restrictionsC;
	delete[] probabilities;

	system((string("notepad ") + wDirectory + string(OUTPUT_F_NAME)).c_str());
}

void findBestLambda(double** phies, size_t** xs, double& lambda, double initial, double step, bool check(double)) {
	double** current_phies = new double* [N];
	size_t* temp = nullptr;
	double xsMax = -INFINITY;

	for (size_t i = 0; i < N; i++) {
		current_phies[i] = new double[Z];
	}

	for (double i = initial; check(i); i += step) {
		double last_max = -INFINITY;

		for (size_t j = 0; j < N; j++) {
			estimatePhies(last_max, last_max, current_phies[j], j, i);
		}

		temp = findXs(current_phies);
		if (validate(temp)) {
			double sum = 0;
			for (size_t i = 0; i < N; i++) {
				sum += temp[i];
			}

			if (sum >= xsMax) {
				*xs = temp;
				xsMax = sum;
				lambda = i;

				for (size_t j = 0; j < N; j++) {
					for (size_t k = 0; k < Z; k++) {
						phies[j][k] = current_phies[j][k];
					}
				}
			}
			else {
				delete[] temp;
			}
		}
		else {
			delete[] temp;
		}
	}
	for (size_t i = 0; i < N; i++) {
		delete[] current_phies[i];
	}
}

void estimatePhies(double& max, const double& last_max, double* phies, size_t x, double lambda) {
	double c_last_max = last_max == -INFINITY ? 0 : last_max;
	double c_max = -INFINITY;

	for (size_t i = 0; i < Z; i++) {
		phies[i] = (1 - pow(probabilities[x], 1 + i)) * exp(lambda * i * restrictionsC[x]) + c_last_max;

		if ((phies[i] < (1 + x)) && phies[i] > c_max && (c_max = phies[i]));
	}

	max = c_max;
};

double estimateResult(size_t* xs) {
	double result = 1.0;

	for (size_t i = 0; i < N; i++) {
		result *= 1 - pow(probabilities[i], 1 + xs[i]);
	}

	return result;
}

size_t* findXs(double** phies) {
	struct {
		size_t x;
		double phi;
	} max_phies[7];

	for (size_t i = 0; i < N; i++) {
		max_phies[i].x = -INFINITY;
		max_phies[i].phi = -INFINITY;

		for (size_t j = 0; j < Z; j++) {
			auto a = phies[i][j];
			if (phies[i][j] > max_phies[i].phi) {
				max_phies[i].x = j;
				max_phies[i].phi = phies[i][j];
			}
		}
	}

	size_t* result = new size_t[N];
	for (size_t i = 0; i < N; i++) {
		result[i] = max_phies[i].x;
	}

	return result;
};

bool validate(size_t* xs) {
	double wS = 0;
	double cS = 0;

	double* temp = new double[7];
	memcpy(temp, xs, 7);

	for (size_t i = 0; i < N; i++) {
		wS += xs[i] * restrictionsW[i];
		cS += xs[i] * restrictionsC[i];
	}

	return wS <= W && cS <= C;
}

void input(ifstream& fs) {
	fs >> N;

	restrictionsW = new double[N];
	restrictionsC = new double[N];
	probabilities = new double[N];

	for (size_t i = 0; i < N; i++) {
		fs >> restrictionsW[i];
	}
	fs >> W;

	for (size_t i = 0; i < N; i++) {
		fs >> restrictionsC[i];
	}
	fs >> C;

	for (size_t i = 0; i < N; i++) {
		fs >> probabilities[i];
		probabilities[i] = 1 - probabilities[i];
	}
}

void printPhies(double** phies, ofstream& fs) {
	for (size_t i = 0; i < Z; i++) {
		string str = to_string(i) += " | ";

		for (size_t j = 0; j < N; j++) {
			str += to_string(phies[j][i]) + ", ";
		}

		str += "\n";

		fs << str.c_str();
	}
};

void showInfo() {
	stringstream infoBuffer;
	infoBuffer << "----------------------------------------------------------------------------------------------\n";
	infoBuffer << "You should have input.txt file where you will write restrictions and items count.\n";
	infoBuffer << "Input format.\n";
	infoBuffer << "\n";
	infoBuffer << "N\n";
	infoBuffer << "c_1 c_2 c_3 ... c_n C\n";
	infoBuffer << "w_1 w_2 w_3 ... w_n W\n";
	infoBuffer << "p_1 p_2 p_3 ... p_n\n";
	infoBuffer << "\n";
	infoBuffer << "Example...\n";
	infoBuffer << "7\n";
	infoBuffer << "4.0 8.0 7.0 3.0 5.0 3.0 6.0 100\n";
	infoBuffer << "4.0 6.0 12.0 10.0 5.0 8.0 10.0 120\n";
	infoBuffer << "0.9 0.85 0.88 0.75 0.9 0.8 0.75\n";
	infoBuffer << "\n";
	infoBuffer << "You can find results in the same directory in output.txt file.\n";
	infoBuffer << "----------------------------------------------------------------------------------------------\n";

	cout << infoBuffer.str();
}
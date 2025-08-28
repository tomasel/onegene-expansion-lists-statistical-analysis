#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <thread>
#include <stack>
#include "metriche.hpp"

int main(int argc, char **argv) {

	// check input
	if (argc != 6) {
		cerr << "This tool needs five arguments: ./a.out <set1> <set2> <dataset1> <dataset2> <savepath>" << endl;
		cerr << "\t<set1>: path to the first set of csv files\n"
						"\t<set2>: path to the second set of csv files\n"
						"\t<dataset1>: csv file containing a list of genes, where for each gene there is an expansion list in <set1>\n"
						"\t<dataset2>: same as <dataset1>, but relative to <set2>"
						"\t<savepath>: where to save the resulting files" << endl;
		return EXIT_FAILURE;
	}
	for (int i = 1; i < argc; ++i) {
		if (!fs::exists(argv[i])) {
			cerr << argv[i] << " does not exist" << endl;
			if (i == 5) {
				if (!fs::create_directories(argv[5])) {
					cerr << "could not create new directory " << argv[5] << endl;
					exit(EXIT_FAILURE);
				}
			} else {
				exit(EXIT_FAILURE);
			}
		}
	}

	fs::path set1 = argv[1];
	fs::path set2 = argv[2];
	fs::path dataset1 = argv[3];
	fs::path dataset2 = argv[4];
	fs::path savepath = argv[5];

	// nome_gene, vettore di file associati a nome_gene
	map<string, fs::path> *dsFiles1 = InizializzaPercorsi(set1);
	map<string, fs::path> *dsFiles2 = InizializzaPercorsi(set2);
	ListaGeni *lista1, *lista2;

	// TEST
	cout << "set1 has " << dsFiles1->size() << " expansions" << endl;
	cout << "set2 has " << dsFiles2->size() << " expansions" << endl;
	int counter_check = 0;

	// parallelizing or not
	const unsigned threadsCount = thread::hardware_concurrency();
	if (threadsCount == 0) {
		cerr << "The number of threads could not be established correctly,"
						"this is going to take a long time!" << endl;
		// do the work
		for (const auto &gene1: *dsFiles1)
			if (auto gene2 = dsFiles2->find(gene1.first); gene2 != dsFiles2->end()) {
				lista1 = LeggiLista(gene1.second);
				lista2 = LeggiLista(gene2->second);
				CalcolaMetricheAll(lista1, lista2, dataset1, dataset2, savepath);
				// free some memory
				delete lista1;
				delete lista2;
				++counter_check;
			}
	} else {
		cout << "using " << threadsCount << " threads" << endl;  // TEST
		// split jobs for every thread
		unsigned curPos = 0;
		vector<vector<vector<fs::path>>> dsFilesParallel(threadsCount);
		/* the idea is this
		 *
		 * dsFilesParallel
		 * └ thread----------thread---------...---thread
		 *   ├ path1, path2  ├ path1, path2       ├...
		 *   ├ path1, path2  ├ path1, path2       ├...
		 *   ...             ...                  ├...
		 *   └ path1, path2  └ path1, path2       └...
		 *
		 * so for every thread we have a vector which can be seen as a list of genes,
		 * and for every gene we have two associated files which will be the input of the function
		 */
		for (auto &it1: *dsFiles1) {
			if (auto it2 = dsFiles2->find(it1.first); it2 != dsFiles2->end()) {
				if (curPos < dsFilesParallel.size()) {
					dsFilesParallel[curPos++].emplace_back(vector<fs::path>{std::move(it1.second), std::move(it2->second)});
				} else {
					curPos = 0;
					dsFilesParallel[curPos++].emplace_back(vector<fs::path>{std::move(it1.second), std::move(it2->second)});
				}
			}
		}
		// spawn threads
		vector<thread> threadsVector;
		for (auto & filesVector : dsFilesParallel) {
			threadsVector.emplace_back(LaunchThreads, dataset1, dataset2, savepath, filesVector);
		}
		for (auto &th: threadsVector) {
				th.join();
		}
	}

	// cleanup
	delete dsFiles1;
	delete dsFiles2;

	return EXIT_SUCCESS;
}

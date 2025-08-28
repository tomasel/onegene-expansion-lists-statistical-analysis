#ifndef METRICHE_H
#define METRICHE_H

#include <ostream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <set>
#include <filesystem>
#include <map>

using namespace std;
namespace fs = std::filesystem;

class Gene {
	string name;
	unsigned rank{};
	float frel{};
	string gene_symb;
	string gene_desc;
public:
	Gene() = default;

	explicit Gene(const string &name);

	Gene(const string &name, unsigned int rank);

	Gene(string name, unsigned int rank, float frel);

	Gene(const string &name, unsigned int rank, float frel, const string &geneSymbol, const string &description);

	const string &getName() const;

	unsigned int getRank() const;

	float getFrel() const;

	void setRank(unsigned int _rank);

	bool operator<(const Gene &) const;

	friend ostream &operator<<(ostream &os, const Gene &gene);
};

class ListaGeni {
public:
	string dataset;                 // dataset di provenienza della lista
	string geneName;                // gene espanso
	string expansionCode;           // codice numerico usato dal database delle espansioni per identificarle univocamente
	string alpha;                   // parametro alpha usato per produrre l'espansione
	set<Gene> *geneSet = nullptr;     // struttura per contenere i dati della lista

	virtual ~ListaGeni();

	friend ostream &operator<<(ostream &os, const ListaGeni &geni);
};

// case-insensitive string comparison
// from https://stackoverflow.com/questions/11635/case-insensitive-string-comparison-in-c
bool ichar_equals(char a, char b);

bool iequals(const std::string &a, const std::string &b);

void LaunchThreads(const fs::path &dataset1,
									 const fs::path &dataset2,
									 const fs::path &savepath,
									 const vector<vector<fs::path>> &filesVector);

vector<string>
tokenizeFirstLineInFile(const fs::path &fp, const char &c);

void OutputHeader(const ListaGeni *lista1,
									const ListaGeni *lista2,
									const fs::path &pathVESP1,
									const fs::path &pathVESP2,
									ostream &os,
									bool inIteration);

ListaGeni *LeggiLista(const fs::path &file);

map<string, fs::path>
*InizializzaPercorsi(const string &path);

void CalcolaMetricheAll(ListaGeni *lista1,
												ListaGeni *lista2,
												const fs::path &pathVESP1,
												const fs::path &pathVESP2,
												const fs::path &pathResults);

int DimensioneIntersezione(const set<Gene> &lista1, const set<Gene> &lista2);

int DimensioneEsclusione(const set<Gene> &lista1, const set<Gene> &lista2);

bool SortGeneRank(const Gene &lista1, const Gene &lista2);

pair<set<Gene>, set<Gene>>
RimuoviNonComuni(const set<Gene> &lista1, const set<Gene> &lista2);

double SpearmanDistance(const set<Gene> &lista1, const set<Gene> &lista2);

double SpearmanRho(const set<Gene> &lista1, const set<Gene> &lista2, bool ties);

double KendallDistance(const set<Gene> &lista1, const set<Gene> &lista2);

double KendallContaGruppiTies(const set<Gene> &lista);

pair<int, int> KendallTies(const set<Gene> &lista1, const set<Gene> &lista2);

int KendallConcordi(const set<Gene> &lista1, const set<Gene> &lista2);

double KendallTauA(const set<Gene> &lista1, const set<Gene> &lista2);

double KendallTauB(const set<Gene> &lista1, const set<Gene> &lista2);

double JaccardDistance(const set<Gene> &lista1, const set<Gene> &lista2);

pair<long double, long double>
FisherTest(const set<Gene> &lista1,
					 const set<Gene> &lista2,
					 const fs::path &pathDataset);

void CalcolaMetriche(const ListaGeni *lista1,
										 const ListaGeni *lista2,
										 const fs::path &pathVESP1,
										 const fs::path &pathVESP2,
										 ostream &os,
										 bool inIteration,
										 float iteration);

void CalcolaCondivisi(const set<Gene> &lista1,
								 const set<Gene> &lista2,
								 const fs::path &geneName,
								 const fs::path &pathListe);

void CalcolaIterativoFreq(const ListaGeni *lista1,
													const ListaGeni *lista2,
													float growth,
													const fs::path &resultPath,
													const fs::path &pathVESP1,
													const fs::path &pathVESP2);

ListaGeni *FirstNFreq(const ListaGeni *lista, float frel);

set<Gene> CorreggiRank(const set<Gene> &lista, bool generaPareggi);

double calculateMean(const set<Gene> &lista);

double calculateStdDev(const set<Gene> &lista);

#endif

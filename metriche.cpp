#include <cassert>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <list>
#include <string>
#include <vector>

#include "fraction.hpp"
#include "metriche.hpp"

using namespace std;
namespace fs = std::filesystem;

// Metodi

Gene::Gene(const string& name,
    unsigned int rank,
    float frel,
    const string& geneSymbol,
    const string& description)
    : name(name)
    , rank(rank)
    , frel(frel)
    , gene_symb(geneSymbol)
    , gene_desc(description)
{
    this->name = name;
    this->rank = rank;
    this->frel = frel;
    this->gene_symb = geneSymbol;
    this->gene_desc = description;
}

bool Gene::operator<(const Gene& gene) const
{
    return name < gene.name;
}

unsigned int Gene::getRank() const
{
    return rank;
}

float Gene::getFrel() const
{
    return frel;
}

Gene::Gene(const string& name, unsigned int rank)
    : name(name)
    , rank(rank)
{
    this->name = name;
    this->rank = rank;
}

Gene::Gene(const string& name)
    : name(name)
{
    this->name = name;
}

Gene::Gene(string name,
    unsigned int rank,
    float frel)
    : name(std::move(name))
    , rank(rank)
    , frel(frel)
{
}

ostream& operator<<(ostream& os, const Gene& gene)
{
    os << gene.name << "," << gene.rank << "," << gene.frel << ","
       << gene.gene_symb << "," << gene.gene_desc;
    return os;
}

void Gene::setRank(unsigned int _rank)
{
    Gene::rank = _rank;
}

const string& Gene::getName() const
{
    return name;
}

ostream& operator<<(ostream& os, const ListaGeni& listaGeni)
{
    os << listaGeni.expansionCode << ","
       << listaGeni.dataset << ","
       << listaGeni.geneName << ","
       << listaGeni.alpha << ","
       << listaGeni.geneSet->size();
    return os;
}

ListaGeni::~ListaGeni()
{
    delete geneSet;
}

// Funzioni helper

bool ichar_equals(char a, char b)
{
    return tolower(static_cast<unsigned char>(a)) == tolower(static_cast<unsigned char>(b));
}

bool iequals(const std::string& a, const std::string& b)
{
    return equal(a.begin(), a.end(), b.begin(), b.end(), ichar_equals);
}

vector<string> tokenizeFirstLineInFile(const fs::path& fp, const char& c)
{
    // stream e strutture di supporto
    vector<string> res;
    string word, temp;
    string dataset, geneName;
    ifstream inStream(fp);

    // legge prima riga e la trasforma in uno stream
    getline(inStream, temp);
    stringstream ssTemp(temp);

    // spezza lo stream in stringhe separate da c
    while (getline(ssTemp, word, c))
        res.push_back(word);

    inStream.close();
    return res;
}

void LaunchThreads(const fs::path& dataset1,
    const fs::path& dataset2,
    const fs::path& savepath,
    const vector<vector<fs::path>>& filesVector)
{
    ListaGeni *lista1, *lista2;
    for (auto& entry : filesVector) {
        lista1 = LeggiLista(entry[0]);
        lista2 = LeggiLista(entry[1]);
        CalcolaMetricheAll(lista1, lista2, dataset1, dataset2, savepath);
        delete lista1;
        delete lista2;
    }
}

ListaGeni* LeggiLista(const fs::path& file)
{
    auto listaGeni = new ListaGeni;
    listaGeni->alpha = "0.05";
    listaGeni->geneSet = new set<Gene>;
    listaGeni->expansionCode = file.filename();
    listaGeni->expansionCode.erase(
        listaGeni->expansionCode.find('_')); // <codice_lista>_<dataset>.<estensione> -> <codice_lista>
    // vettore di stringhe per contenere i singoli valori letti dal file .csv
    vector<string> tmp = tokenizeFirstLineInFile(file, ' ');
    listaGeni->dataset = tmp[0];
    listaGeni->geneName = tmp[3];
    for (auto& c : listaGeni->geneName) { // "VIT_..." -> "vit_..."
        c = (char)tolower(c);
    }
    // correggo il nome del gene, potrebbe avere il prefisso "vv_"
    // se la lista è indicata come Vv_23 va indicato come la prima espansione con "_1"
    if (listaGeni->geneName.rfind("vv_", 3) != string::npos) {
        if (iequals(listaGeni->dataset, "Vv_23")) {
            listaGeni->dataset.append("_1");
        }
        listaGeni->geneName.erase(0, 3);
    }

    ifstream sLista(file);
    string temp, word;
    // salto le prime due righe (nomi delle colonne)
    for (int i = 0; i < 2; ++i)
        getline(sLista, temp);
    // leggo i dati
    while (getline(sLista, temp)) {
        tmp.clear();
        stringstream ssTemp(temp);
        while (getline(ssTemp, word, ',')) {
            tmp.push_back(word);
        }
        listaGeni->geneSet->insert(Gene(tmp[1], stoi(tmp[0]), stof(tmp[3])));
    }
    sLista.close();
    return listaGeni;
}

map<string, fs::path>* InizializzaPercorsi(const string& path)
{
    auto* dsFiles = new map<string, fs::path>();
    string geneName;
    vector<string> tmp;
    for (const auto& file : fs::directory_iterator(path)) {
        tmp = tokenizeFirstLineInFile(file, ' ');
        assert(tmp.size() == 4);
        geneName = tmp[3];
        for (auto& c : geneName) // "VIT_..." -> "vit_..."
            c = (char)tolower(c);
        if (geneName.rfind("vv_", 3) != string::npos) // se inizia con "vv_"; da C++20 è possibile usare starts_with()
            geneName.erase(0, 3);
        (*dsFiles)[geneName] = file.path();
    }
    return dsFiles;
}

void CalcolaMetricheAll(ListaGeni* lista1,
    ListaGeni* lista2,
    const fs::path& pathVESP1,
    const fs::path& pathVESP2,
    const fs::path& pathResults)
{
    // directory per i risultati
    fs::path gene_dir = pathResults / lista1->geneName;
    if (!fs::create_directories(gene_dir)) {
        cerr << "could not create " << gene_dir << endl;
        exit(EXIT_FAILURE);
    }
    gene_dir /= lista1->expansionCode + "_" + lista2->expansionCode + "-frel.csv";
    if (fs::exists(gene_dir)) {
        cerr << gene_dir.filename() << " already exists, skipping..." << endl;
        return;
    }
    // calcola le metriche iterando sul valore frel (decrescente)
    CalcolaIterativoFreq(lista1, lista2, 0.01, gene_dir, pathVESP1, pathVESP2);
}

int DimensioneIntersezione(const set<Gene>& lista1, const set<Gene>& lista2)
{
    int n = 0;
    for (const auto& gene : lista1)
        if (lista2.find(gene) != lista2.end())
            ++n;
    return n;
}

int DimensioneEsclusione(const set<Gene>& lista1, const set<Gene>& lista2)
{
    int n = (int)lista1.size();
    for (auto& gene : lista2)
        if (lista1.find(gene) != lista1.end())
            --n;
    return n;
}

bool SortGeneRank(const Gene& lista1, const Gene& lista2)
{
    return lista1.getRank() < lista2.getRank();
}

pair<set<Gene>, set<Gene>> RimuoviNonComuni(const set<Gene>& lista1, const set<Gene>& lista2)
{
    set<Gene> m1 = lista1;
    set<Gene> m2 = lista2;
    for (auto& gene : lista1)
        if (auto tmp = lista2.find(gene); tmp == lista2.end()) {
            m1.erase(gene);
        }
    for (auto& gene : lista2)
        if (auto tmp = lista1.find(gene); tmp == lista1.end()) {
            m2.erase(gene);
        }
    assert(m1.size() == m2.size());
    return { m1, m2 };
}

double calculateMean(const set<Gene>& lista)
{
    double sum = 0;
    for (const auto& gene : lista) {
        sum += gene.getRank();
    }
    return sum / (double)lista.size();
}

double calculateStdDev(const set<Gene>& lista)
{
    double mean = calculateMean(lista);
    double sumSquaredDiff = 0;
    for (const auto& gene : lista) {
        sumSquaredDiff += pow(gene.getRank() - mean, 2);
    }
    return dsqrtl(sumSquaredDiff / (double)lista.size());
}

double calculateCovariance(const set<Gene>& lista1, const set<Gene>& lista2, size_t size)
{
    double sum = 0;
    double Ex;
    double Ey;
    auto x_m = calculateMean(lista1);
    auto y_m = calculateMean(lista2);
    for (auto& x : lista1) {
        auto y = lista2.find(x);
        Ex = (double)x.getRank() - x_m;
        Ey = (double)y->getRank() - y_m;
        sum += Ex * Ey;
    }
    return sum / (double)size;
}

int KendallConcordi(const set<Gene>& lista1, const set<Gene>& lista2)
{
    int concordi = 0;
    for (auto i = lista1.begin(); i != lista1.end(); ++i) {
        for (auto j = next(i); j != lista1.end(); ++j) {
            auto x = lista2.find(*i);
            auto y = lista2.find(*j);
            if (x != lista2.end() && y != lista2.end()) {
                // se la coppia è concorde
                if ((i->getRank() > j->getRank() && x->getRank() > y->getRank()) || (i->getRank() < j->getRank() && x->getRank() < y->getRank())) {
                    ++concordi;
                }
            }
        }
    }
    return concordi;
}

double KendallContaGruppiTies(const set<Gene>& lista)
{
    double ties = 1; // numero di elementi a parimerito nell'i-esimo gruppo di pareggi (sempre almeno uno)
    double sumTies = 0;
    for (auto it = lista.begin(); it != lista.end(); ++it) {
        if (next(it) != lista.end()) {
            if (it->getRank() == next(it)->getRank()) {
                ++ties;
            } else {
                sumTies += ties * (ties - 1) / 2;
                ties = 1;
            }
        } else {
            sumTies += ties * (ties - 1) / 2;
        }
    }
    return sumTies;
}

pair<int, int> KendallTies(const set<Gene>& lista1, const set<Gene>& lista2)
{
    // this function assumes lista1 and lista2 to be ordered by ranks
    unsigned T = 0; // counter for tied elements with equal rank
    unsigned ties1 = 0; // ties in lista1
    // check lista1
    for (auto it = lista1.begin(); it != lista1.end(); ++it) {
        if (auto nextit = next(it); (it->getRank() == nextit->getRank()) && nextit != lista1.end()) {
            ++T;
        } else {
            if (T >= 1)
                ties1 += coeff_binom(T + 1, 2).compute2(); // T did not count the last value, so we add 1
            T = 0;
        }
    }
    // check lista2
    unsigned ties2 = 0; // ties in lista2
    T = 0; // ensure counter is reset
    for (auto it = lista2.begin(); it != lista2.end(); ++it) {
        if (auto nextit = next(it); (it->getRank() == nextit->getRank()) && nextit != lista2.end()) {
            ++T;
        } else {
            if (T >= 1) {
                ties2 += coeff_binom(T + 1, 2).compute2();
                T = 0;
            }
        }
    }
    return make_pair(ties1, ties2);
}

ListaGeni* FirstNFreq(const ListaGeni* lista, float frel)
{
    auto* res = new ListaGeni(*lista);
    res->geneSet = new set<Gene>;
    for (auto& g : *lista->geneSet)
        if (g.getFrel() >= frel)
            res->geneSet->insert(g);
    return res;
}

void OutputHeader(const ListaGeni* lista1,
    const ListaGeni* lista2,
    const fs::path& pathVESP1,
    const fs::path& pathVESP2,
    ostream& os,
    bool inIteration)
{
    if (inIteration)
        os << "Frel,";
    else
        os << "Lista,";
    os << "Geni condivisi,"
          "Dimensione lista "
       << lista1->expansionCode << "," << "Dimensione lista " << lista2->expansionCode << "," << "Spearman distance,"
                                                                                                 "Kendall distance,"
                                                                                                 "Jaccard distance,"
                                                                                                 "Spearman rho (no ties),"
                                                                                                 "Spearman rho (ties),"
                                                                                                 "Kendall tau-a (no ties),"
                                                                                                 "Kendall tau-b (ties),"
                                                                                                 "Fisher test ("
       << (string)pathVESP1.filename() << "),pvalue," << "Fisher test (" << (string)pathVESP2.filename() << "),pvalue" << endl;
}

set<Gene> CorreggiRank(const set<Gene>& lista, bool generaPareggi)
{
    std::list<Gene> tmp;
    set<Gene> res;
    for (auto& gene : lista)
        tmp.push_back(gene);
    tmp.sort(SortGeneRank);
    int rank = 1;
    for (auto& gene : tmp) { // rank contigui mantenendo l'ordine
        gene.setRank(rank++);
    }
    if (generaPareggi) {
        int group = 1; // gruppi con rank uguale, si parte dal gruppo con rank = 1
        for (auto it = tmp.begin(); it != tmp.end(); ++it) {
            auto nextit = std::next(it);
            if (it->getFrel() - nextit->getFrel() >= 0.010 && nextit != tmp.end()) {
                it->setRank(group);
                ++group;
            } else {
                it->setRank(group);
            }
        }
        for (auto gene : tmp) {
            res.insert(std::move(gene));
        }
    } else {
        for (auto gene : tmp) {
            res.insert(std::move(gene));
        }
    }
    return res;
}

// Funzioni per il calcolo delle metriche

double SpearmanDistance(const set<Gene>& lista1, const set<Gene>& lista2)
{
    int distance = 0;
    if (lista1.empty() || lista2.empty() || DimensioneIntersezione(lista1, lista2) == 0)
        return NAN;
    for (auto& gene : lista1)
        if (auto tmp = lista2.find(gene); tmp != lista2.end())
            distance += abs((int)(gene.getRank() - tmp->getRank()));
    return distance;
}

double SpearmanRho(const set<Gene>& lista1, const set<Gene>& lista2, bool ties)
{
    auto n = lista1.size();
    if (n < 2)
        return NAN;
    if (ties) {
        double num = calculateCovariance(lista1, lista2, n);
        double den = calculateStdDev(lista1) * calculateStdDev(lista2);
        return num / den;
    }
    unsigned sqdistance = 0;
    for (auto& gene1 : lista1)
        if (auto gene2 = lista2.find(gene1); gene2 != lista2.end())
            sqdistance += (unsigned)pow(gene1.getRank() - gene2->getRank(), 2);
    fraction res;
    res.mul(6);
    res.mul(sqdistance);
    res.div(n);
    res.div((unsigned)(pow(n, 2) - 1));
    return (double)(1.0 - res.compute2());
}

double JaccardDistance(const set<Gene>& lista1, const set<Gene>& lista2)
{
    unsigned intersection = DimensioneIntersezione(lista1, lista2);
    if (lista1.empty() || lista2.empty() || intersection == 0)
        return 1; // max distance if not comparable
    unsigned join = lista1.size();
    for (auto& gene : lista2)
        if (auto tmp = lista1.find(gene); tmp == lista1.end())
            ++join;
    return (intersection == join) ? 0 : 1 - intersection / (double)join;
}

pair<long double, long double> FisherTest(const set<Gene>& lista1,
    const set<Gene>& lista2,
    const fs::path& pathDataset)
{
    // values to return
    fraction ftest;
    long double pvalue = 0;
    set<Gene> setDataSet;
    set<Gene> setUnione;
    ifstream inDataset(pathDataset);
    // initialize gene lists
    for (auto& gene : lista1) {
        setUnione.insert(Gene(gene.getName()));
    }
    for (auto& gene : lista2) {
        setUnione.insert(Gene(gene.getName()));
    }
    // initialize dataset
    vector<string> row;
    string word, temp;
    string dataset, geneName;
    getline(inDataset, temp); // skip first line
    while (getline(inDataset, temp)) {
        temp.erase(temp.find(',')); // gene name is the first value, discard the rest of the line
        for (auto& c : temp) // "VIT_..." -> "vit_..."
            c = (char)tolower(c);
        setDataSet.insert(Gene(temp));
    }
    inDataset.close();
    // initialize table values
    int aAndb = DimensioneIntersezione(lista1, lista2);
    int aNotb = DimensioneEsclusione(lista1, lista2);
    int bNota = DimensioneEsclusione(lista2, lista1);
    int allNotunion = DimensioneEsclusione(setDataSet, setUnione);
    // p-value
    int marginalColumn1 = aAndb + aNotb;
    int marginalRow1 = aAndb + bNota;
    int N = (int)setDataSet.size();
    long double observedProb = hg1(aAndb, marginalRow1, marginalColumn1, N).compute2();
    long double updatedf;
    // sum the probabilites for all the possible tables with fixed parameters
    for (int j = min(marginalRow1, marginalColumn1); j >= (aAndb > allNotunion ? allNotunion : 0); --j) {
        updatedf = hg1(j, marginalRow1, marginalColumn1, N).compute2();
        if (updatedf <= observedProb) {
            pvalue += updatedf;
        }
    }
    // fisher test
    if (aNotb == 0 || bNota == 0) {
        /* Adding 0.5 to all cells is commonly done to improve the small sample properties
         * from the Chi-square test (and it asymptotically removes the first-order bias from the estimate
         * of the log-odds ratio)
         * https://stats.stackexchange.com/questions/198571/dealing-with-0-in-cell-count-for-fishers-exact-test*/
        // num = (aAndb + 0.5) / (aNotb + 0.5);
        // den = (bNota + 0.5) / (allNotunion + 0.5);
        ftest.mul(2 * aAndb + 1);
        ftest.div(2 * aNotb + 1);
        ftest.mul(2 * allNotunion + 1);
        ftest.div(2 * bNota + 1);
    } else {
        // num = (double)aAndb / aNotb;
        // den = (double)bNota / allNotunion;
        ftest.mul(aAndb);
        ftest.div(aNotb);
        ftest.mul(allNotunion);
        ftest.div(bNota);
    }
    return make_pair(ftest.compute2(), pvalue);
}

double KendallDistance(const set<Gene>& lista1, const set<Gene>& lista2)
{
    int distance = 0;
    for (auto it11 = lista1.begin(); it11 != lista1.end(); ++it11) {
        for (auto it12 = next(it11); it12 != lista1.end(); ++it12) {
            auto it21 = lista2.find(*it11);
            auto it22 = lista2.find(*it12);
            if (it21 != lista2.end() && it22 != lista2.end()) {
                // se la coppia è "non concorde" aumento la distanza di 1
                if ((it11->getRank() > it12->getRank() && it21->getRank() < it22->getRank()) || (it11->getRank() < it12->getRank() && it21->getRank() > it22->getRank())) {
                    ++distance;
                }
            }
        }
    }
    return distance;
}

double KendallTauA(const set<Gene>& lista1, const set<Gene>& lista2)
{
    auto n = DimensioneIntersezione(lista1, lista2);
    if (lista1.empty() || lista2.empty() || n == 0)
        return NAN;
    auto nd = KendallDistance(lista1, lista2);
    auto n0 = (int)coeff_binom(n, 2).compute2();
    auto nc = KendallConcordi(lista1, lista2);
    auto res = (double)(nc - nd) / n0;
    return res;
}

// Maurice G. Kendall, “The treatment of ties in ranking problems”, Biometrika Vol. 33, No. 3, pp. 239-251. 1945
double KendallTauB(const set<Gene>& lista1, const set<Gene>& lista2)
{
    if (lista1.empty() || lista2.empty() || DimensioneIntersezione(lista1, lista2) == 0)
        return NAN;
    auto Q = KendallDistance(lista1, lista2); // coppie discordi
    auto P = KendallConcordi(lista1, lista2); // coppie concordi
    auto pairTies = KendallTies(lista1, lista2);
    auto T = pairTies.first;
    auto U = pairTies.second;
    auto n = coeff_binom(lista1.size(), 2).compute2();
    return (double)(P - Q) / (dsqrtl((n - T) * (n - U))); // utilizza la definizione corretta
}

// Wrapper

void CalcolaMetriche(const ListaGeni* lista1,
    const ListaGeni* lista2,
    const fs::path& pathVESP1,
    const fs::path& pathVESP2,
    ostream& os, bool inIteration,
    float iteration)
{
    // save some useful values for later
    unsigned dimInt = DimensioneIntersezione(*lista1->geneSet, *lista2->geneSet);
    string geneName = lista1->geneName;
    unsigned l1size = lista1->geneSet->size();
    unsigned l2size = lista2->geneSet->size();
    // bool comparable = !(lista1->geneSet->empty() || lista2->geneSet->empty() || dimInt == 0);

    // remove non-shared genes and correct ranks without generating ties
    auto returnPair = RimuoviNonComuni(*lista1->geneSet, *lista2->geneSet);
    auto l1noties = CorreggiRank(returnPair.first, false);
    auto l2noties = CorreggiRank(returnPair.second, false);

    // generate lists with ties
    auto l1ties = CorreggiRank(returnPair.first, true);
    auto l2ties = CorreggiRank(returnPair.second, true);

    // calculate metrics
    auto Sd = (unsigned)SpearmanDistance(l1noties, l2noties);
    auto Kd = (unsigned)KendallDistance(l1noties, l2noties);
    auto Sr_noties = SpearmanRho(l1noties, l2noties, false);
    auto Sr_ties = SpearmanRho(l1ties, l2ties, true);
    auto Ka = KendallTauA(l1noties, l2noties);
    auto Kb = KendallTauB(l1ties, l2ties);
    pair<long double, long double> fisher1 = FisherTest(*lista1->geneSet, *lista2->geneSet, pathVESP1);
    pair<long double, long double> fisher2 = FisherTest(*lista1->geneSet, *lista2->geneSet, pathVESP2);
    auto Jd = JaccardDistance(*lista1->geneSet, *lista2->geneSet);

    // print everything to file
    if (inIteration) {
        os << iteration << ",";
    } else {
        OutputHeader(lista1, lista2, pathVESP1, pathVESP2, os, false);
        os << geneName << ",";
    }
    os << dimInt << "," << l1size << "," << l2size << "," << Sd << "," << Kd << "," << Jd << "," << Sr_noties << "," << Sr_ties << "," << Ka << "," << Kb << "," << fisher1.first << "," << fisher1.second << "," << fisher2.first << "," << fisher2.second << endl;
}

void CalcolaCondivisi(const set<Gene>& lista1,
    const set<Gene>& lista2,
    const fs::path& geneName,
    const fs::path& pathListe)
{
    // per recuperare le descrizioni dei geni è necessario estrarle dalle liste di espansione
    // quindi leggo tutte le liste e aggiungo la descrizione
    //  - controllo stringa in row[1]
    //  - se hit, carico e associo row[9]
    string outfile_name;
    outfile_name.append(geneName);
    outfile_name.erase(outfile_name.find('.'));
    outfile_name.append("-condivisi.csv");
    ofstream os(outfile_name, ostream::out | ostream::app);
    vector<string> row;
    string word, temp;
    vector<string> files;
    set<Gene> genesFound; // ottimizzazione: struttura d'appoggio per controllare quali geni sono già stati letti
    set<Gene> intersection;

    os << "Nome gene,rank 2021,rank 2023,frel 2021,frel 2023,Gene symbol,Description,Delta Frel" << endl;

    // creo l'insieme intersezione e ricalcolo i rank
    for (auto& gene : lista1) {
        if (lista2.find(gene) != lista2.end()) {
            intersection.insert(gene);
        }
    }
    // creo una lista di tutti i file
    for (const auto& file : filesystem::directory_iterator(pathListe))
        files.push_back(file.path());
    // cerco gene symbol e gene_desc nelle liste di adiacenza
    for (const auto& file : files) {
        // se ho già trovato tutti i geni, ho finito
        if (genesFound.size() < intersection.size()) {
            ifstream infile(file, iostream::in);
            // salto le prime due righe di ogni file
            // contengono i nomi delle colonne e vanno scartati
            for (int i = 0; i < 2; ++i) {
                std::getline(infile, temp);
            }
            while (getline(infile, temp)) {
                row.clear();
                for (auto& c : temp)
                    c = (char)tolower(c);
                stringstream ss(temp);
                while (getline(ss, word, ',')) {
                    row.push_back(word);
                }
                const auto& found = intersection.find(Gene(row[1]));
                // se trovo il gene nell'intersezione E se non l'ho mai letto prima
                if (found != intersection.end() && genesFound.find(Gene(row[1])) == genesFound.end()) {
                    auto it1 = lista1.find(Gene(row[1]));
                    auto it2 = lista2.find(Gene(row[1]));
                    assert(it1 != lista1.end() && it2 != lista2.end());
                    genesFound.insert(Gene(row[1])); // segnalo il gene trovato
                    std::string name = found->getName();
                    unsigned rank = it1->getRank();
                    unsigned rank2 = it2->getRank();
                    float frel = it1->getFrel();
                    float frel2 = it2->getFrel();
                    // accedere ad un elemento out of bound in un vector ha comportamento non definito
                    std::string symbol = (row.size() < 5) ? "NaN" : row[4];
                    std::string desc = (row.size() < 10) ? "NaN" : row[9];
                    float delta_frel = frel - frel2;
                    os << name << ","
                       << rank << ","
                       << rank2 << ","
                       << frel << ","
                       << frel2 << ","
                       << symbol << ","
                       << desc << ","
                       << delta_frel << endl;
                }
            }
        }
    }
}

void CalcolaIterativoFreq(const ListaGeni* lista1,
    const ListaGeni* lista2,
    float growth,
    const fs::path& resultPath,
    const fs::path& pathVESP1,
    const fs::path& pathVESP2)
{
    float iteration = 1; // frel di partenza
    ofstream outFile(resultPath, ostream::out | ostream::app);
    OutputHeader(lista1, lista2, pathVESP1, pathVESP2, outFile, true);
    while (iteration >= 0) {
        CalcolaMetriche(FirstNFreq(lista1, iteration),
            FirstNFreq(lista2, iteration),
            pathVESP1,
            pathVESP2,
            outFile,
            true,
            iteration);
        iteration = roundf(iteration * 100) / 100 - growth;
    }
}

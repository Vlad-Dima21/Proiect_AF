#include <bits/stdc++.h>

using namespace std;

class Triplet   //folositor la APM
{
public:
    int a;  //parintele nodudlui din heap in APM
    int b;  //nodul curent din heap
    int c;  //costul cu care s-a ajuns la nodul curent
};


class CompareAPM  //clasa pentru priority queue din APM
{
public:
    bool operator() (Triplet a, Triplet b)
    {
        return a.c > b.c;   // compar costurile
    }
};

class ComparePair // -----------//------------- din Dijkstra
{
public:
    bool operator() (pair<int,int> a, pair<int,int> b)
    {
        return a.second > b.second;   // compar costurile
    }
};

class Graf
{
    int nrNoduri;
    int nrMuchii;

    vector<vector<int>> lsAd; ///lista de adiacenta a grafului
    vector<vector<pair<int,int>>> lsAdCost; // lista de adiacenta a grafului ponderat

    ///metode private
    void CompBicTimp(const int nodCrt, int &nrComp, int* parinte,int* timpIntr, int* timpMin, stack<pair<int,int>>* stiva, set<int>* noduriComp);///metoda care calculeaza timpii fiecarui nod

    void CTCCalc(const int nodCrt, int &nrComp,int* timpIntr, int* timpMin, stack<int>* stiva, vector<int>* noduriComp, bool* inStiva);///metoda care calculeaza efectiv timpii pentru CTC

    void SortTopDFS(const int nodCrt, list<int> &listaTop, bool* vizitat);

    void CritConDFS(const int nodCrt,int* timpIntr,int* timpMin,int* parinte,vector<vector<int>>& muchiiCrit);

    int DisjointInterogareCompresie(const int x, vector<int>& tata);
    void DisjointReuniune(const int x, const int y, vector<int>& tata, vector<int>& rang);

public:
    Graf():nrNoduri(0),nrMuchii(0) {}///constructor pentru problemele care nu au date de intrare in format standard(initializare prin metoda)
    Graf(bool orientat,ifstream &in);
    ~Graf();

    void BFS(char* fisierIntrare ="bfs.in",char* fisierIesire={"bfs.out"});///primul exercitiu(afiseaza intr-un fisier lungimea drumurilor de la un nod la restul)

    void DFS(char* fisierIesire="dfs.out");

    void CompBic(ofstream &out);///metoda principala a exercitiului componente biconexe(face dfs din noduri)


    void CTC(ofstream &out);///metoda principala a exercitiului CTC(Componente tare conexe)


    void SortTop(ofstream &out);///metoda pentru sortarea topologica


    void HavelHakimi();

    void CritConCreateGraph(int n, vector<vector<int>>& connections, vector<vector<int>>& muchiiCrit);///metoda pentru leetcode CCN(cu tot cu crearea listei de adiacenta)
    
    void APM(ifstream &in, vector<vector<int>>& muchiiApm, int &costMin);///metoda pentru APM

    void Disjoint(ifstream& in, ofstream& out, const int N, const int M, vector<int>& tata, vector<int>& rang);

    void Dijkstra(ifstream& in, vector<int>& drum);

    bool BellmanFord(ifstream& in, vector<int>& cost); //intoarce true daca s-a putut aplica algoritmul
                                                       //adica daca nu contine un ciclu, altfel false



};

Graf::~Graf()
{
    lsAd.clear();
}

Graf::Graf(bool orientat,ifstream &in)
{
    int x, y;

    in>>nrNoduri>>nrMuchii;
    lsAd.resize(nrNoduri+1);

    for (int i = 0; i < nrMuchii; i++)
        {
            in>>x>>y;
            lsAd[x].push_back(y);
            if (!orientat)
               lsAd[y].push_back(x);
        }
}
void Graf::DFS(char* fisierIesire)
{
    stack<int> stiva;

    int nrComponente = 0;
    bool* nodMarcat = new bool[nrNoduri+1]();///initializat cu false pentru toate nodurile

    for (int i = 1; i <= nrNoduri; i++)
    {
        if (!nodMarcat[i])
            {
                stiva.push(i);
                nodMarcat[i] = true;
                int nodCrt;

                while (!stiva.empty())
                {
                    nodCrt = stiva.top();
                    stiva.pop();

                    for (auto nodAd : lsAd[nodCrt])
                    {
                        if (!nodMarcat[nodAd])
                        {
                            nodMarcat[nodAd] = true;
                            stiva.push(nodAd);
                        }
                    }
                }
                nrComponente++;
            }
    }
    delete[] nodMarcat;
    ///scrierea in fisier
    ofstream out(fisierIesire);
    out<<nrComponente;
    out.close();
}


void Graf::BFS(char* fisierIntrare, char* fisierIesire)
{
    int s_bfs;  ///nodul de pornire
    ifstream in(fisierIntrare);
    int x, y;

    in>>nrNoduri>>nrMuchii>>s_bfs;
    lsAd.resize(nrNoduri+1); ///un fel de initializare a listei de adiacenta in constructor

    for (int i = 0; i < nrMuchii; i++)
    {
        in>>x>>y;
        lsAd[x].push_back(y);
    }
    in.close();

    queue<int> coada;

    bool* nodVizitat = new bool[nrNoduri+1];///array-ul cu care verific daca a fost vizitat un nod
    int* distNod = new int[nrNoduri+1];///nr muchii catre fiecare nod

    for (int i = 1; i <= nrNoduri; i++)
        {
            distNod[i] = (i == s_bfs ? 0 : -1);
            nodVizitat[i] = false;
        }

    coada.push(s_bfs);
    int nodCrt;///nodul curent din parcurgerea laterala
    while(!coada.empty())
    {
        nodCrt = coada.front();
        coada.pop();
        nodVizitat[nodCrt] = true;

        for (int i = 0; i < lsAd[nodCrt].size(); i++)///trec prin toate nodurile adiacente
        {
            if (!nodVizitat[lsAd[nodCrt][i]])///daca nodul este nevizitat il adaug in coada
            {
                distNod[lsAd[nodCrt][i]] = distNod[nodCrt]+1;
                nodVizitat[lsAd[nodCrt][i]] = true;
                coada.push(lsAd[nodCrt][i]);
            }
        }
    }

    ///scrierea in fisier
    ofstream out(fisierIesire);

    for (int i = 1; i <= nrNoduri; i++)
        out<<distNod[i]<<' ';
    out.close();
    delete[] distNod;
    delete[] nodVizitat;
}

void Graf::CompBic(ofstream &out)
{
    int nrComp = 0; ///variabila in care memorez numarul de componente biconexe gasite
    set<int>* noduriComp = new set<int>[nrNoduri];///am mai multe multimi in care sunt nodurile fiecarei componente

    int* parinte = new int[nrNoduri+1]();///array cu parintii fiecarui nod(pentru a verifica daca o muchie este intre tata si un fiu)

    int* timpIntr = new int[nrNoduri+1]();///timpii de intrare in noduri
    int* timpMin = new int[nrNoduri+1]();///timpul minim la care se poate ajunge(eventual printr-o muchie de intoarcere)
    stack<pair<int,int>>* stiva = new stack<pair<int,int>>;///stiva folosita

    ///array-urile dinamice sunt deja initializate cu 0 pe toate pozitiile

    for (int i = 1; i <= nrNoduri; i++)
    {
        if (timpIntr[i] == 0)///nu am vizitat nodul inca
            {
                CompBicTimp(i,nrComp,parinte,timpIntr,timpMin,stiva,noduriComp);
                if (!stiva->empty()) ///dupa apelul functiei au ramas noduri(componenta nu contine muchie de intoarcere)
                {

                    nrComp++;
                    while (!stiva->empty())
                    {
                        pair<int,int> m = stiva->top();
                        stiva->pop();
                        noduriComp[nrComp-1].insert(m.first);
                        noduriComp[nrComp-1].insert(m.second);
                    }
                }
            }
    }
    ///partea de scriere in fisier
    out<<nrComp<<'\n';

    for (int i = 0; i < nrComp; i++)
        {
            for (auto itr : noduriComp[i])
                out<<itr<<' ';
            out<<'\n';
        }
    delete[] noduriComp;
    delete[] parinte;
    delete[] timpIntr;
    delete[] timpMin;
    delete stiva;
}

void Graf::CompBicTimp(const int nodCrt, int &nrComp, int* parinte,int* timpIntr, int* timpMin, stack<pair<int,int>>* stiva, set<int>* noduriComp)
{
    static int timp = 1;///variabila statica pentru calcularea timpilor de intrare

    timpIntr[nodCrt] = timpMin[nodCrt] = timp++;

    for (auto fiu : lsAd[nodCrt])
    {
        if (timpIntr[fiu] == 0)
        {
            parinte[fiu] = nodCrt;
            stiva->push(make_pair(nodCrt,fiu));///adaugam muchia in stiva

            CompBicTimp(fiu,nrComp,parinte,timpIntr,timpMin,stiva,noduriComp);

            timpMin[nodCrt] = min(timpMin[nodCrt], timpMin[fiu]); ///actualizez timpul(nivelul) minim la care poate ajunge nodul curent

            if (timpMin[fiu] >= timpIntr[nodCrt])///fiul nu poate ajunge la un stramos al lui nodCrt(deci avem o componenta)
            {
                nrComp++;

                pair<int,int> p;
                while (!(stiva->top().first == nodCrt && stiva->top().second == fiu))
                {
                    p = stiva->top();
                    stiva->pop();
                    noduriComp[nrComp-1].insert(p.first);
                    noduriComp[nrComp-1].insert(p.second);
                }
                p = stiva->top();
                stiva->pop();
                noduriComp[nrComp-1].insert(p.first);
                noduriComp[nrComp-1].insert(p.second);
            }
        }
        else if (fiu != parinte[nodCrt])///am gasit o muchie de intoarcere
            timpMin[nodCrt] = min(timpMin[nodCrt], timpIntr[fiu]);
    }
}

void Graf::CTC(ofstream &out)
{
    ///asemanator cu biconex doar ca ma intereseaza doar nodurile si nu muchiile

    int nrComp = 0;
    int* parinte = new int[nrNoduri+1]();
    int* timpIntr = new int[nrNoduri+1]();
    int* timpMin = new int[nrNoduri+1]();
    stack<int>* stiva = new stack<int>;
    vector<int>* noduriComp = new vector<int>[nrNoduri];///incerc cu vector in loc de set ca nu mai am nevoie si e si mai rapid
    bool* inStiva = new bool[nrNoduri+1]();///in plus fata de biconex, intr-o componenta tare conexa pot fi doar nodurile din dfs

    for (int i = 1; i <= nrNoduri; i++)
        if (timpIntr[i] == 0)///daca nodul este nevizitat incep dfs
            CTCCalc(i,nrComp,timpIntr,timpMin,stiva,noduriComp,inStiva);

    out<<nrComp<<'\n';
    for (int i = 0; i < nrComp; i++)
    {
        for (auto itr : noduriComp[i])
            out<<itr<<' ';
        out<<'\n';
    }

    delete[] inStiva;
    delete[] noduriComp;
    delete stiva;
    delete[] timpMin;
    delete[] timpIntr;
    delete[] parinte;
}

void Graf::CTCCalc(const int nodCrt, int &nrComp,int* timpIntr, int* timpMin, stack<int>* stiva, vector<int>* noduriComp, bool* inStiva)
{
    static int timp = 1;

    timpIntr[nodCrt] = timpMin[nodCrt] = timp++;

    stiva->push(nodCrt);
    inStiva[nodCrt] = true;

    for (auto fiu : lsAd[nodCrt])
    {
        if (timpIntr[fiu] == 0) ///fiu nevizitat
        {
            inStiva[fiu] = true;

            CTCCalc(fiu,nrComp,timpIntr,timpMin,stiva,noduriComp,inStiva);

            timpMin[nodCrt] = min(timpMin[nodCrt], timpMin[fiu]);
        }
        else if (inStiva[fiu])///fiul este in stiva(muchia este de intoarcere)
            timpMin[nodCrt] = min(timpMin[nodCrt], timpIntr[fiu]);
    }
    if (timpIntr[nodCrt] == timpMin[nodCrt])///nodul curent este radacina componentei tare conexe
        {
            nrComp++;

            while (stiva->top() != nodCrt)///toate nodurile din stiva pana la nodul radacina fac parte din ctc
            {
                inStiva[stiva->top()] = false;
                noduriComp[nrComp-1].push_back(stiva->top());
                stiva->pop();
            }
            inStiva[stiva->top()] = false;
            noduriComp[nrComp-1].push_back(stiva->top());
            stiva->pop();
        }
}

void Graf::SortTop(ofstream &out)
{
    list<int> listaTop;///lista co nodurile sortate topologic

    bool* vizitat = new bool[nrNoduri+1]();

    for (int i = 1; i <= nrNoduri; i++)
    {
        if (vizitat[i] == false)
            SortTopDFS(i,listaTop,vizitat);
    }

    ///scriere fisier

    for (auto itr : listaTop)
        out<<itr<<' ';

    delete[] vizitat;
}

void Graf::SortTopDFS(const int nodCrt, list<int> &listaTop, bool* vizitat)
{
    vizitat[nodCrt] = true;
    for (auto fiu : lsAd[nodCrt])
        if (!vizitat[fiu])
            {
                SortTopDFS(fiu,listaTop,vizitat);
                vizitat[fiu] = true;///pentru fiii comuni
            }

    ///adaug nodul curent atunci cand se iese din el in timpul parcurgerii dfs
    listaTop.push_front(nodCrt);
}


void Graf :: HavelHakimi()
    {
        int nr_noduri;

        cout<<"Algoritm HavelHakimi\n\n";
        cout<<"Numarul de noduri: ";
        cin>>nr_noduri;

        int* grade = new int[nr_noduri];

        cout<<"Grade noduri: ";

        for (int i = 0; i < nr_noduri; i++)
            cin>>grade[i];

        for (int i = 0; i < nr_noduri; i++)
        {
            sort(grade+i, grade + nr_noduri, greater<int>());///sortam descrescator gradele

//            ///debugging
//            for (int j = i; j < nr_noduri; j++)
//                cout<<grade[j]<<' ';
//            cout<<'\n';
//            ///debugging

            if (nr_noduri-i-1 < 0)///nodul are gradul mai mare decat numarul de noduri ramase(imposibil)
                {
                    cout<<"Gradele nu pot fi ale unui graf simplu";
                    return;
                }

            int contor_zero = 0;///pentru a verifica daca toate elementele au gradul 0(ceea ce inseamna ca formeaza un graf)

            for (int j = i+1; j <= i + grade[i]; j++)
            {

                if (grade[j] == 0)///nodul i are cel putin o muchie incidenta careia ii lipsete al doilea nod
                    {
                        cout<<"Gradele nu pot fi ale unui graf simplu";
                        return;
                    }
                grade[j]--;
                if (grade[j] == 0)
                    contor_zero++;
            }
            if (contor_zero == nr_noduri-i-1)///restul gradelor sunt 0 deci exista graful
            {
                cout<<"Gradele pot fi ale unui graf simplu";
                return;
            }
        }
        delete[] grade;
    }

void Graf :: CritConCreateGraph(const int n, vector<vector<int>>& connections, vector<vector<int>>& muchiiCrit)
{
    lsAd.resize(n);
    nrMuchii = n;

    for (auto muchie : connections)
        {
            lsAd[muchie[0]].push_back(muchie[1]);
            lsAd[muchie[1]].push_back(muchie[0]);
        }

    ///principiu asemanator cu biconex si ctc

    int* timpIntr = new int[n]();
    int* timpMin = new int[n]();
    int* parinte = new int[n]();

    for (int i = 0; i < n; i++)
    {
        if (timpIntr[i] == 0)
            CritConDFS(i,timpIntr,timpMin,parinte,muchiiCrit);
    }

    delete[] parinte;
    delete[] timpMin;
    delete[] timpIntr;
}

void Graf :: CritConDFS(const int nodCrt,int* timpIntr,int* timpMin,int* parinte,vector<vector<int>>& muchiiCrit)
{
    static int timp = 1;

    timpIntr[nodCrt] = timpMin[nodCrt] = timp++;

    for (auto fiu : lsAd[nodCrt])
    {
        if (timpIntr[fiu] == 0)
        {
            parinte[fiu] = nodCrt;
            CritConDFS(fiu,timpIntr,timpMin,parinte,muchiiCrit);

            timpMin[nodCrt] = min(timpMin[nodCrt], timpMin[fiu]);

            if (timpMin[fiu] > timpIntr[nodCrt])///daca fiul nu poate ajunge la un stramos al tatalui inseamna ca muchia dintre cei doi e critica
            {
                muchiiCrit.push_back({nodCrt,fiu});
            }
        }
        else if(parinte[nodCrt] != fiu)
            timpMin[nodCrt] = min(timpMin[nodCrt], timpIntr[fiu]);
    }
}

void Graf :: APM(ifstream &in, vector<vector<int>>& muchiiApm, int &costMin)
{
    // Partea de initializare a grafului
    in>>nrNoduri>>nrMuchii;

    lsAdCost.resize(nrNoduri+1);

    int x, y, c;

    for (int i = 0; i < nrMuchii; i++)
        {
            in>>x>>y>>c;
            lsAdCost[x].push_back(make_pair(y,c));
            lsAdCost[y].push_back(make_pair(x,c));
        }
    
    //Algoritmul propriu-zis
    costMin = 0;

    priority_queue<Triplet,vector<Triplet>,CompareAPM> heap; //folosesc priority queue drept heap

    vector<bool> vizitat(nrNoduri+1); //daca e vizitat inseamna ca este deja in apm

    heap.push(Triplet{0,1,0}); 

    while (!heap.empty())
    {
        Triplet curent = heap.top();
        heap.pop();

        // in triplet:
        //  a -> parintele nodului curent
        //  b -> nodul curent la care s-a ajuns
        //  c -> costul cu care s-a ajuns la nodul curent
        if (!vizitat[curent.b])
        {
            vizitat[curent.b] = true;

            costMin += curent.c;

            if (curent.b != 1)
            {
                muchiiApm.push_back({curent.a,curent.b});
            }

            for (auto vecin : lsAdCost[curent.b])
                {
                    heap.push(Triplet{curent.b,vecin.first,vecin.second});
                }
        }
    }
}

void Graf :: Disjoint(ifstream& in, ofstream& out, const int N, const int M, vector<int>& tata, vector<int>& rang)
{
    int cod, x, y;

    for (int i = 0; i < M; i++)
        {
            in >> cod >> x >> y;

            if (cod == 1)   //reunirea multimilor
                DisjointReuniune(DisjointInterogareCompresie(x,tata),DisjointInterogareCompresie(y,tata),tata,rang); //metoda care face reuniunea
            else
                {
                    if (DisjointInterogareCompresie(x,tata) == DisjointInterogareCompresie(y,tata)) //metoda returneaza radacina arborelui multime(facand si compresia drumurilor)
                        out<<"DA"<<'\n';
                    else
                        out<<"NU"<<'\n';
                }
        }
}

int Graf :: DisjointInterogareCompresie(const int x, vector<int>& tata)
{
    int radacina = x; //aflu radacina arborelui

    while (radacina != tata[radacina]) //conditia de oprire e ca tatal nodului sa fie el insusi(adica cand se ajunge la radacina)
        radacina = tata[radacina];

    int nodCrt = x; //nodul curent din drumul de la x la radacina

    while (nodCrt != tata[nodCrt]) //incep sa fac compresia drumului
        {
            int aux = tata[nodCrt];
            tata[nodCrt] = radacina;
            nodCrt = aux;
        }

    return radacina;
}

void Graf :: DisjointReuniune(const int x, const int y, vector<int>& tata, vector<int>& rang)
{
    if (rang[x] > rang[y])
        {
            tata[y] = x;
            rang[y] = rang[x];
        }
    else if (rang[x] == rang[y])
        {
            tata[x] = y;
            rang[x] = ++rang[y];
        }
    else
        {
            tata[x] = y;
            rang[x] = rang[y];
        }
}

void Graf :: Dijkstra(ifstream& in, vector<int>& drum)
{
    //partea de initializare
    in >> nrNoduri >> nrMuchii;
    lsAdCost.resize(nrNoduri+1);
    drum.resize(nrNoduri+1);

    for (int i = 2; i <= nrNoduri; i++) //nodurile care nu sunt accesibile din 1 au lungimea 0 (asa e in cerinta)
        drum[i] = 0;

    int x, y, cost;

    for (int i = 0; i < nrMuchii; i++)
    {
        in >> x >> y >> cost;
        lsAdCost[x].push_back(make_pair(y,cost));
    }

    //algoritm propriu-zis

    priority_queue<pair<int,int>,vector<pair<int,int>>,ComparePair> heap; //folosesc din nou un priority queue drept heap

    vector<bool> vizitat(nrNoduri+1);

    heap.push(make_pair(1,0));

    while (!heap.empty())
    {
        pair<int,int> tupluCrt = heap.top(); //in tuplu am nodul si costul catre nod
        heap.pop();

        if (!vizitat[tupluCrt.first])
        {
            vizitat[tupluCrt.first] = true;
            drum[tupluCrt.first] = tupluCrt.second;

            for (auto vecin : lsAdCost[tupluCrt.first])
            {
                //pentru fiecare vecin il adaug in heap cu costul curent plus costul arcului
                heap.push(make_pair(vecin.first,tupluCrt.second + vecin.second));
            }
        }
    }
}

bool Graf :: BellmanFord(ifstream& in, vector<int>& cost)
{
    in >> nrNoduri >> nrMuchii;
    lsAdCost.resize(nrNoduri + 1);
    cost.resize(nrNoduri + 1);

    int x, y, c;
    for (int i = 0; i < nrMuchii; i++)
    {
        in >> x >> y >> c;
        lsAdCost[x].push_back(make_pair(y,c));
    }
    
    queue<int> coada;
    vector<int> nrImbunatatiri(nrNoduri+1); //pastrez pentru fiecare nod de cate ori i-a fost imbunatatit costul
                                            //daca aceasta valoare a depasit nrNoduri - 1 inseamna ca in graf exista
                                            //un ciclu negativ   
    
    vector<bool> inCoada(nrNoduri+1);       //in coada adaug doar nodurile al caror cost s-a putut imbunatati

    bool areCicluNeg = false;               //primeste true cand este indeplinita conditia de mai sus

    for (int i = 1; i <= nrNoduri; i++)
    {
        cost[i] = 1000000000;
    }
    cost[1] = 0;

    coada.push(1);
    inCoada[1] = true;

    while (!coada.empty() && !areCicluNeg)
    {
        int nodCrt = coada.front();
        coada.pop();

        inCoada[nodCrt] = false;
        for (auto vecin : lsAdCost[nodCrt])
        {
            if (cost[vecin.first] > cost[nodCrt] + vecin.second)    //vecin.first -> nodul catre care duce arcul // vecin.second -> costul arcului
                if (nrImbunatatiri[vecin.first] + 1 > nrNoduri) //daca este satisfacuta conditia de mai sus...
                {
                    areCicluNeg = true;
                    break;
                }
                else
                {
                    cost[vecin.first] = cost[nodCrt] + vecin.second;
                    coada.push(vecin.first);
                    inCoada[vecin.first] = true;
                    nrImbunatatiri[vecin.first] += 1;
                }
        }
    }

    return !areCicluNeg;
}

///pentru leetcode
class Solution {
public:
    vector<vector<int>> criticalConnections(int n, vector<vector<int>>& connections) {
        vector<vector<int>> muchiiCrit;

        Graf g;
        g.CritConCreateGraph(n,connections,muchiiCrit);

        return muchiiCrit;
    }
};



int main()
{
    ///pentru BFS

//    Graf g;
//    g.BFS();

    ///pentru DFS

//    ifstream in("dfs.in");
//    Graf g(false,in);
//    g.DFS("dfs.out");
//    in.close();

    ///pentru Componente Biconexe

//    ifstream in("biconex.in");
//    ofstream out("biconex.out");
//    Graf g(false,in);
//    in.close();
//    g.CompBic(out);
//    out.close();

    ///pentru Componente Tare Conexe

//    ifstream in("ctc.in");
//    Graf g(true, in);
//    ofstream out("ctc.out");
//    g.CTC(out);
//    in.close();
//    out.close();

    ///pentru Sortare Topologica

//    ifstream in("sortaret.in");
//    Graf g(true,in);
//    ofstream out("sortaret.out");
//    g.SortTop(out);
//    in.close();
//    out.close();

    ///pentru Havel Hakimi
    
//    Graf g;
//    g.HavelHakimi();

    // pentru APM

// ifstream in("apm.in");
// ofstream out("apm.out");
// vector<vector<int>> muchiiApm;
// int costMin;
// Graf g;
// g.APM(in,muchiiApm,costMin); //returneaza prin parametrii
// in.close();
// out<<costMin<<'\n';
// out<<muchiiApm.size()<<'\n';
// for (auto muchie : muchiiApm)
//     out<<muchie[0]<<' '<<muchie[1]<<'\n';
// out.close();


// pentru Disjoint

{
    // ifstream in("disjoint.in");
    // ofstream out("disjoint.out");

    // int N, M;
    
    // in >> N >> M;

    // vector<int> tata(N+1), rang(N+1);

    // for (int i = 1; i <= N; i++)
    //     {
    //         tata[i] = i;
    //         rang[i] = 1;
    //     }
    
    // Graf g;
    // g.Disjoint(in,out,N,M,tata,rang);   //citirea operatiilor si afisarea se realizeaza in metoda
    
    // in.close();
    // out.close();


}

// pentru Dijkstra

{
    // ifstream in("dijkstra.in");
    // ofstream out("dijkstra.out");

    // vector<int> drum; //in el 'returnez' rezultatul

    // Graf g;
    // g.Dijkstra(in,drum);
    // in.close();

    // for (int i = 2; i < drum.size(); i++)
    // {
    //     out << drum[i] << ' ';
    // }
    // out.close();
}

// pentru Bellman-Ford
{
    ifstream in("bellmanford.in");
    ofstream out("bellmanford.out");

    vector<int> cost;

    Graf g;
    bool merge = g.BellmanFord(in,cost);

    if (merge)
    {
        for (int i = 2; i < cost.size(); i++)
            out << cost[i] << ' ';
    }
    else
    {
        out << "Ciclu negativ!";
    }
}

}

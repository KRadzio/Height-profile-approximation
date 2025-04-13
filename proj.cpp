// g++ proj.cpp -o proj
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

#define ENDL '\n'
#define TEXT ".txt"
#define NODES_NAME "nodesValues"
#define EQ 4
#define POINTS_NUM 512
#define LAGRANGE 0
#define SPLINES 1

struct point
{
    double x;
    double y;
};

// Obsługa plików
vector<point> loadFile(string fileName, bool skipFirstLine = true, char separator = ',')
{
    ifstream file(fileName);
    string line;
    vector<point> v;
    if (skipFirstLine)
        getline(file, line);
    while (getline(file, line))
    {
        point p;
        stringstream ss(line);
        string element;
        getline(ss, element, separator);
        p.x = stod(element);
        getline(ss, element, ENDL);
        p.y = stod(element);
        v.push_back(p);
    }
    return v;
}
// używane przy faktoryzacji LU
void swapRows(double **U, double **L, double **permMatrix, int index, int N)
{
    double pivot = abs(U[index][index]);
    int pivotRow = index;

    for (int i = index + 1; i < N; i++)
    {
        if (abs(U[i][index]) > pivot)
        {
            pivot = abs(U[i][index]);
            pivotRow = i;
        }
    }
    // zamiana rzędów
    if (pivotRow != index)
    {
        double tmp;
        for (int i = 0; i < N; i++)
        {
            if (i >= index)
            {
                tmp = U[index][i];
                U[index][i] = U[pivotRow][i];
                U[pivotRow][i] = tmp;
            }
            else
            {
                tmp = L[index][i];
                L[index][i] = L[pivotRow][i];
                L[pivotRow][i] = tmp;
            }
            tmp = permMatrix[index][i];
            permMatrix[index][i] = permMatrix[pivotRow][i];
            permMatrix[pivotRow][i] = tmp;
        }
    }
}
// faktoryzacja LU zadziała zawsze i da dokładne wyniki(dane w tym zadaniu nie są zbyt duże więc można jej użyć)
void LU(double **A, double *b, double *res, int N)
{
    // setup
    double **U = new double *[N];          // górna trójkątna część
    double **L = new double *[N];          // dolna trójkątna część
    double **permMatrix = new double *[N]; // do rozwiązania układu
    double *y = new double[N];             // L * y = b w podstawieniu w przód

    for (int i = 0; i < N; i++)
    {
        U[i] = new double[N];
        L[i] = new double[N];
        permMatrix[i] = new double[N];
        y[i] = 0;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                L[i][j] = permMatrix[i][j] = 1;
            U[i][j] = A[i][j];
        }
    }
    for (int i = 0; i < N - 1; i++) // kolumny
    {
        swapRows(U, L, permMatrix, i, N);
        for (int j = i + 1; j < N; j++) // wiersze
        {
            L[j][i] = U[j][i] / U[i][i];
            for (int k = i; k < N; k++) // k zastępuje numer kolumny
                U[j][k] = U[j][k] - L[j][i] * U[i][k];
        }
    }

    for (int i = 0; i < N; i++) // podstawienie w przód
    {
        double element = 0.0;

        for (int j = 0; j < i; j++)
            element += L[i][j] * y[j];

        y[i] = (b[i] - element) / L[i][i];
    }

    for (int i = N - 1; i >= 0; i--) // podstawienie w tył
    {
        double element = 0.0;

        for (int j = i + 1; j < N; j++)
            element += U[i][j] * res[j]; // res jest wypełniony 1

        res[i] = (y[i] - element) / U[i][i];
    }
    // res zawiera wektor rozwiązań
    delete y;
    for (int i = 0; i < N; i++)
    {
        delete U[i];
        delete L[i];
        delete permMatrix[i];
    }
    delete U;
    delete L;
    delete permMatrix;
}

// metody interpolacji
double Lagrange(vector<point> &points, double currDist)
{
    double yApp = 0.0;
    for (int i = 0; i < points.size(); i++)
    {
        double currElement = 1.0;
        for (int j = 0; j < points.size(); j++)
            if (i != j)
                currElement *= (currDist - points[j].x) / (points[i].x - points[j].x);

        yApp += currElement * points[i].y;
    }
    return yApp;
}

double splines(vector<point> &points, double currDist)
{
    int nodesNum = points.size();
    int eqNum = EQ * (nodesNum - 1); // n węzłów, n-1 układów równań
    // są w "zwykłych" tablicach a nie w wektorach ze względu na to że w poprzdnim projekcie faktoryzacja LU była zaimplemantowana na takich tablicach
    double **EqMatrix = new double *[eqNum]; // macierz zawierająca wszyskie n-1 układów równań
    double *coeffsVec = new double[eqNum];   // a0, b0, c0, d0, ....... an-1, bn-1, cn-1, dn-1
    double *res = new double[eqNum];         // wyniki mnożenia EqMatrix * coeffsVec

    for (int i = 0; i < eqNum; i++)
    {
        EqMatrix[i] = new double[eqNum];
        coeffsVec[i] = 0;
        res[i] = 1; // monożenie przez 1 nie powoduje zmian
    }

    for (int i = 0; i < eqNum; i++) // wypełnienie macierzy zerami
        for (int j = 0; j < eqNum; j++)
            EqMatrix[i][j] = 0;

    // S0(x) = a0 + b0(x − x0) + c0(x − x0)^2 + d0(x − x0)^3
    // S1(x) = a1 + b1(x − x1) + c1(x − x1)^2 + d1(x − x1)^3
    // S0′(x) = b0 + 2c0(x − x0) + 3d0(x − x0)^2
    // S1′(x) = b1 + 2c1(x − x1) + 3d1(x − x1)^2
    // S0′′(x) = 2c0 + 6d0(x − x0)
    // S1′′(x) = 2c1 + 6d1(x − x1)

    // uzupełnienie pierwszego równania
    double h = points[1].x - points[0].x;

    // (1) S0 = f(x0) ; a0 = f(x0) ; (x0 - x0) = 0 -> co zeruje pozostałe współyczniki w tym przypadku
    EqMatrix[0][0] = 1;
    coeffsVec[0] = points[0].y;

    // (2) S0(x1) = f(x1); a0 + b0*h + c0*h^2 + d0*h^3
    EqMatrix[1][0] = 1;
    EqMatrix[1][1] = h;
    EqMatrix[1][2] = pow(h, 2);
    EqMatrix[1][3] = pow(h, 3);
    coeffsVec[1] = points[1].y;

    // (7) i (8) S''0(x0) = 0, bo S0 = a0
    EqMatrix[2][2] = 1;
    coeffsVec[2] = 0;

    // (6) S''n-1(xn) = 0
    h = points[nodesNum - 1].x - points[nodesNum - 2].x;
    EqMatrix[3][4 * (nodesNum - 2) + 2] = 2;
    EqMatrix[3][4 * (nodesNum - 2) + 3] = 6 * h;
    coeffsVec[3] = 0;

    for (int i = 1; i < nodesNum - 1; i++)
    {
        h = points[i].x - points[i - 1].x;

        // (3) Si(xi) = f(xi) ; ai = f(xi)
        EqMatrix[4 * i][4 * i] = 1;
        coeffsVec[4 * i] = points[i].y; // ai

        // (4) Si(xi+1) = f(x+1) ; ai + bi*h + ci*h^2 + di*h^3 = f(xi+1)
        EqMatrix[4 * i + 1][4 * i] = 1;
        EqMatrix[4 * i + 1][4 * i + 1] = h;
        EqMatrix[4 * i + 1][4 * i + 2] = pow(h, 2);
        EqMatrix[4 * i + 1][4 * i + 3] = pow(h, 3);
        coeffsVec[4 * i + 1] = points[i + 1].y; // bi

        // (5)/* S'i-1(xi) = S'i(xi) ; bi-1 + 2ci-1*h + 3di-1*h^2  - bi = 0
        EqMatrix[4 * i + 2][4 * (i - 1) + 1] = 1;
        EqMatrix[4 * i + 2][4 * (i - 1) + 2] = 2 * h;
        EqMatrix[4 * i + 2][4 * (i - 1) + 3] = 3 * pow(h, 2);
        EqMatrix[4 * i + 2][4 * i + 1] = -1;
        coeffsVec[4 * i + 2] = 0; // ci

        // (6) S''i-1(xi) = S''i(xi) ; 2ci-1 + 6di-1*h - 2ci = 0
        EqMatrix[4 * i + 3][4 * (i - 1) + 2] = 2;
        EqMatrix[4 * i + 3][4 * (i - 1) + 3] = 6 * h;
        EqMatrix[4 * i + 3][4 * i + 2] = -2;
        coeffsVec[4 * i + 3] = 0; // di
    }
    LU(EqMatrix, coeffsVec, res, eqNum);

    double yApp;
    for (int i = 0; i < nodesNum - 1; i++)
    {
        yApp = 0;
        if (currDist >= points[i].x && currDist <= points[i + 1].x)
        {
            for (int j = 0; j < 4; j++)
            {
                h = currDist - points[i].x;
                yApp += res[4 * i + j] * pow(h, j);
            }
            break;
        }
    }

    delete res;
    delete coeffsVec;
    for (int i = 0; i < eqNum; i++)
        delete EqMatrix[i];
    delete EqMatrix;

    return yApp;
}

void interpolate(vector<point> &points, vector<point> &interpolationNodes, string filePath)
{
    string methodL = "Lagrange";
    string methodS = "Splines";

    ofstream file;
    ofstream file2;
    file.open(filePath + methodS + TEXT);
    file2.open(filePath + methodL + TEXT);

    double valueApp;
    int stepSize = points[1].x; // przesuwamy się o "średnią" odległość pomiędzy punktami pomiarowymi, currDist jest odległością od początku trasy
    for (int currDist = interpolationNodes[0].x; currDist <= interpolationNodes[interpolationNodes.size() - 1].x; currDist += stepSize)
    {
        valueApp = splines(interpolationNodes, currDist); // można by raz zrobić macierz bo wywołanie różnią się tylko wartością "currDist" (chyba)
        file << currDist << " " << valueApp << endl;
        valueApp = Lagrange(interpolationNodes, currDist);
        file2 << currDist << " " << valueApp << endl;
    }
    file.close();
    file2.close();

    file.open(filePath + NODES_NAME + TEXT); // te same wartości w węzłach dla obu metod

    for (int i = 0; i < interpolationNodes.size(); i++)
        file << interpolationNodes[i].x << " " << interpolationNodes[i].y << endl;

    file.close();
}

void setUpNodesVector(vector<point> &nodes, vector<point> &data, int nodesNum)
{
    int diff = POINTS_NUM / (nodesNum - 1); // 5, 9, 17, 33, 65 -> 4, 8, 16, 32, 64; aby policzyć wartości dla całego przedziału
    int pos = 0;
    for (int j = 0; j < nodes.size(); j++)
    {
        nodes[j].x = data[pos].x;
        nodes[j].y = data[pos].y;
        pos += diff;
        if (pos == POINTS_NUM)
            pos--;
    }
}

int main()
{
    vector<point> mountEverest = loadFile("MountEverest.csv");
    vector<point> tczewStarogard = loadFile("tczew_starogard.txt", false, ' ');
    vector<point> ulmLugano = loadFile("ulm_lugano.txt", false, ' ');
    vector<point> gdansk = loadFile("SpacerniakGdansk.csv");
    vector<string> fileNames = {"MountEverest", "TczewStarogard", "UlmLugano", "SpacerniakGdansk"};
    string outputDirName = "./output/";
    // węzłów jest tak napradę n + 1
    vector<int> nodesNumbers = {8, 16, 32, 64, 128}; // każda próbka zawiera 512 punktów więc używamy potęg dwójki
    vector<point> interpolationNodes;

    for (int i = 0; i < nodesNumbers.size(); i++)
    {
        int nodesNum = POINTS_NUM / nodesNumbers[i] + 1; // inaczej nie będziemy mieli kawałka przedziału
        interpolationNodes = vector<point>(nodesNum);

        setUpNodesVector(interpolationNodes, mountEverest, nodesNum);
        interpolate(mountEverest, interpolationNodes, outputDirName + fileNames[0] + to_string(nodesNum));

        setUpNodesVector(interpolationNodes, tczewStarogard, nodesNum);
        interpolate(tczewStarogard, interpolationNodes, outputDirName + fileNames[1] + to_string(nodesNum));

        setUpNodesVector(interpolationNodes, ulmLugano, nodesNum);
        interpolate(ulmLugano, interpolationNodes, outputDirName + fileNames[2] + to_string(nodesNum));

        setUpNodesVector(interpolationNodes, gdansk, nodesNum);
        interpolate(gdansk, interpolationNodes, outputDirName + fileNames[3] + to_string(nodesNum));

        interpolationNodes.clear();
    }

    return 0;
}